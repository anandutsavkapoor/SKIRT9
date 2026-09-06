/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASLINEEMISSION_HPP
#define GASLINEEMISSION_HPP

#include "Basics.hpp"
#include "StoredTable.hpp"
#include "StoredTableDictionary.hpp"
#include <functional>
class SimulationItem;

//////////////////////////////////////////////////////////////////////

/** GasLineEmission provides species-agnostic building blocks for gas line emission, shared by the
    material mixes:

    - Line inventory: 26 built-in lines (H and He recombination lines, optical forbidden metal
      lines) indexed by LineIndex, plus the extended inventory appended at setup by
      extendLineRegistry().
    - Recombination lines (Case B): L = eps(T, n_e) n_e n_ion V from emissivity tables (Storey &
      Sochi 2015 for H I, Porter et al. 2012 for He I, Storey & Hummer 1995 for He II) loaded by
      initializeRecombinationTables(); without tables the H lines use the legacy P_B form (Storey
      & Hummer 1995; Hui & Gnedin 1997; McClymont, Smith & Tacchella 2025) and the He lines are
      zero. Lyman-alpha always uses the legacy form (the Case B tables exclude the Lyman series).
    - Collisional lines in the nebular limit (level populations set by electron collisions at the
      local T and n_e, no radiative pumping): the statistical-equilibrium solver on atomic models
      loaded by initializeAtomicModels(), or the legacy precomputed q_col(T, n_e) tables when the
      atomic data are unavailable.
    - Statistical equilibrium: the general level-population solver (collisions with any set of
      partners, spontaneous decay, optional radiative pumping) on atomic models read from the
      NonLTELineGasMix species files. The solver follows Kosei Matsumoto's NonLTELineGasMix
      implementation (Matsumoto et al. 2023) and is shared by NonLTELineGasMix and
      DiffuseIonizedGasMix.

    Functions that are pure computations on their arguments (hydrogenLineLuminosity(),
    loadAtomicModel(), solveLevelPopulations(), lineEmissivities()) remain static and need no
    instance. Everything else reads or writes one of this class's three registries (the line
    registry itself, the Case B emissivity-table registry, and the atomic-model registry) and is
    therefore instance state: each consumer that calls initializeRecombinationTables(),
    initializeAtomicModels() or extendLineRegistry() owns its own GasLineEmission instance. This
    also means these registries, and the resources they hold open (the Case B
    StoredTableDictionary and the StoredTable instances spawned from it), are released
    deterministically when that owning instance is destroyed, rather than at unspecified static
    destruction order. */
class GasLineEmission final
{
public:
    /** The constructor initializes the line registry with the built-in lines (the first numLines
        entries, in LineIndex order); the extended inventory is empty until extendLineRegistry()
        is called. */
    GasLineEmission();

    /** The destructor releases any resources acquired by initializeRecombinationTables() (the
        Case B StoredTableDictionary and the StoredTable instances spawned from it). */
    ~GasLineEmission();

    /** The copy constructor and copy assignment operator are deleted because an instance may own
        open resources (see initializeRecombinationTables()) that cannot be meaningfully shared by
        naive member-wise copying. */
    GasLineEmission(const GasLineEmission&) = delete;
    GasLineEmission& operator=(const GasLineEmission&) = delete;

    // ============== Line inventory ==============

    // Line indices for the emission line array
    enum LineIndex {
        // H recombination lines (Case B)
        Lya = 0,   // Lyman-alpha   1215.67 A
        Ha,        // Balmer-alpha   6562.80 A
        Hb,        // Balmer-beta    4861.33 A
        Hg,        // Balmer-gamma   4340.46 A
        Hd,        // Balmer-delta   4101.73 A
        HeBalmer,  // Balmer-epsilon 3970.07 A
        Paa,       // Paschen-alpha  18751 A
        Pab,       // Paschen-beta   12818 A
        Bra,       // Brackett-alpha 40512 A

        // He recombination lines (Case B)
        HeI5876,   // He I  5876 A
        HeI6678,   // He I  6678 A
        HeI7065,   // He I  7065 A
        HeI10830,  // He I  10830 A
        HeII1640,  // He II 1640 A
        HeII4686,  // He II 4686 A

        // Collisional lines (optical forbidden metal lines)
        NII6548,   // [NII] 6548
        NII6583,   // [NII] 6583
        OI6300,    // [OI] 6300
        OI6364,    // [OI] 6364
        OII3729,   // [OII] 3729  (2D 5/2)
        OII3726,   // [OII] 3726  (2D 3/2)
        OIII4363,  // [OIII] 4363
        OIII4959,  // [OIII] 4959
        OIII5007,  // [OIII] 5007
        SII6716,   // [SII] 6716
        SII6731,   // [SII] 6731
    };
    static constexpr int numLines = 26;

    /** One registry line. carrierIonIndex/elementIndex address the consumer's ion-fraction and
        abundance arrays (PhotoIonizationSolver layout); both are -1 for recombination lines, whose
        recombining ion (H II, He II, He III) is selected by the family flags. */
    struct LineDef
    {
        double wavelength;    // rest-frame wavelength [m]
        double mass;          // particle mass for Doppler broadening [kg]
        int carrierIonIndex;  // PhotoIonizationSolver ion stage index, -1 for recombination lines
        int elementIndex;     // abundance element index (0=C..7=Fe), -1 for recombination lines
        bool isHeIRecomb;     // He I recombination line (recombining ion He II)
        bool isHeIIRecomb;    // He II recombination line (recombining ion He III)
    };

    /** The line registry: the first numLines entries are the built-in lines (LineIndex order),
        followed by any lines added by extendLineRegistry(). */
    const std::vector<LineDef>& lineRegistry() const;

    /** One atomic species for the extended inventory: resource base name (e.g. "N_II") plus the
        consumer's carrier indices. */
    struct SpeciesSpec
    {
        std::string name;     // atomic data file base name
        int carrierIonIndex;  // ion stage index in the consumer's ion fraction array
        int elementIndex;     // abundance element index
    };

    /** Appends every recombination table listed in the per-set wavelength index files and every
        radiative transition of each loadable species (all resolved via FilePaths::resource();
        absent species are skipped; built-in lines are not duplicated). Requires
        initializeRecombinationTables() and initializeAtomicModels() to have completed (throws
        FatalError otherwise). Idempotent. The \em item argument is passed through to
        loadAtomicModel() for each newly loaded species (used only to retrieve the units system and
        a logger); this function itself issues a single log message summarizing how many new atomic
        models were loaded, regardless of how many species were considered. Returns the number of
        lines added. */
    int extendLineRegistry(const SimulationItem* item, const std::vector<SpeciesSpec>& species);

    // The rest-frame wavelength, mass, and (for collisional lines) carrier ion/element index of
    // each built-in line are pure implementation detail of the constructor and a few internal
    // helpers; see lineWavelengths, lineMasses, lineCarrierIonIndex, and lineElementIndex near the
    // top of the .cpp. lineRegistry() is the public way to get this information.

    // ============== Recombination lines ==============

    /** Legacy H line luminosity [W] (Lya through Bra): h nu P_B(T, n_e) gammaHI nHI V, with P_B
        the Case B probability that a recombination emits the line. Densities in cm^-3, volume in
        cm^3. */
    static double hydrogenLineLuminosity(int lineIdx, double T, double ne, double gammaHI, double nHI, double V_cm3);

    /** H or He recombination line luminosity [W] (Lya through HeII4686): eps(T, ne) ne nIon V from
        the emissivity tables when loaded, else the legacy form for H (which uses gammaHI, nHI) and
        zero for He. nIon is the recombining ion density (H II, He II or He III). Densities in
        cm^-3, volume in cm^3. */
    double recombinationLineLuminosity(int lineIdx, double T, double ne, double nIon, double gammaHI, double nHI,
                                       double V_cm3) const;

    /** Loads the Case B emissivity tables (StoredTable, resolved via FilePaths::resource()) and
        keeps the resulting StoredTable instances open for as long as this instance exists (they
        are memory-mapped, so this holds no data in memory beyond what is actually looked up). The
        item argument is passed through to StoredTable::open() for each table (used only to
        retrieve a logger; not otherwise a SimulationItem in this hierarchy). Idempotent; throws
        FatalError on a missing or malformed file. */
    void initializeRecombinationTables(const SimulationItem* item);

    /** Returns true when initializeRecombinationTables() has completed. */
    bool recombinationTablesReady() const;

    // ============== Collisionally excited lines ==============

    /** Nebular-limit line luminosity [W] for a collisional line (NII6548 through SII6731, or an
        extended-inventory line): level populations from electron collisions
        at (T, ne) without radiative pumping, times nIon V. Uses the atomic models when loaded,
        else the legacy q_col tables. Densities in cm^-3, volume in cm^3. */
    double collisionalLineLuminosity(int lineIdx, double T, double ne, double nIon, double V_cm3) const;

    /** Loads the atomic models of the built-in collisional lines' carrier species and maps each
        line to its transition by nearest wavelength. Idempotent. The \em item argument is passed
        through to loadAtomicModel() for each species (used only to retrieve the units system and a
        logger). */
    void initializeAtomicModels(const SimulationItem* item);

    /** Returns true when initializeAtomicModels() has completed. */
    bool atomicModelsReady() const;

    /** Model slot of the species carrying the given registry line, or -1 if the line is not served
        by a loaded atomic model. Lines with the same slot share one solve. */
    int lineModelSlot(int lineIdx) const;

    /** Transition index of the given registry line within its model (valid when lineModelSlot() >= 0). */
    int lineTransition(int lineIdx) const;

    // ============== Statistical equilibrium ==============

    /** Static atomic data for one species in SI units, as read by loadAtomicModel(). */
    struct AtomicModel
    {
        double mass{0.};                  // particle mass [kg]
        std::vector<double> energy;       // level energies [J]
        std::vector<double> weight;       // statistical weights
        std::vector<int> indexUpRad;      // upper level index per radiative transition
        std::vector<int> indexLowRad;     // lower level index per radiative transition
        std::vector<double> einsteinA;    // Einstein A coefficient [1/s]
        std::vector<double> einsteinBul;  // Einstein B_ul (per-wavelength convention)
        std::vector<double> einsteinBlu;  // Einstein B_lu (per-wavelength convention)
        std::vector<double> center;       // line center wavelength [m]
        std::vector<double> branchRatio;  // A divided by the sum of A from the same upper level

        /** Collisional transition data for one collision partner. */
        struct ColPartner
        {
            std::string name;                      // partner name ("e-", "H", ...)
            std::vector<double> T;                 // temperature grid [K]
            std::vector<int> indexUpCol;           // upper level index per collisional transition
            std::vector<int> indexLowCol;          // lower level index per collisional transition
            std::vector<std::vector<double>> Kul;  // de-excitation rate [m3/s] on the T grid
        };
        std::vector<ColPartner> colPartner;

        int numLevels() const { return static_cast<int>(energy.size()); }
        int numLines() const { return static_cast<int>(einsteinA.size()); }
    };

    /** Per-cell inputs for solveLevelPopulations(). An empty meanJ skips radiative pumping. The
        optional warn callback receives solver warnings; the solver never logs. */
    struct Environment
    {
        double Tkin{0.};                               // kinetic temperature [K]
        double nTotal{0.};                             // species number density [m^-3]
        std::vector<double> nPartner;                  // collision partner densities [m^-3]
        std::vector<double> meanJ;                     // per-line mean intensity; empty -> pumping skipped
        std::function<void(const std::string&)> warn;  // optional warning sink
    };

    /** Fills model from the resource files NAME_Mass.txt, NAME_Energy.txt, NAME_Rad_Coeff.txt and
        NAME_Col_PARTNER_Temp.txt / NAME_Col_PARTNER_Coeff.txt (NAME the species name, PARTNER each
        collision partner), each resolved via TextInFile as a resource file, keeping at most
        maxNumLevels levels. Throws FatalError if a file cannot be found. The \em item argument is
        passed through to TextInFile for each of these files (used only to retrieve the units
        system and a logger); every file is opened silently (see TextInFile), so a caller loading
        many species in a loop should issue its own single summary log message for the operation as
        a whole. */
    static void loadAtomicModel(const SimulationItem* item, const std::string& speciesName,
                                const std::vector<std::string>& partnerNames, int maxNumLevels, AtomicModel& model);

    /** The atomic model loaded by initializeAtomicModels() into the given slot (see lineModelSlot()). */
    const AtomicModel& atomicModel(int slot) const;

    /** Solves the statistical-equilibrium rate matrix and returns level populations [m^-3]
        normalized to env.nTotal. Throws FatalError if the matrix is singular or the
        solution is not finite. */
    static std::vector<double> solveLevelPopulations(const AtomicModel& model, const Environment& env);

    /** Line power densities [W m^-3] for all radiative transitions: n_up A h c / lambda. */
    static std::vector<double> lineEmissivities(const AtomicModel& model, const std::vector<double>& pops);

    // ================== Data members ==================

private:
    // the mutable line registry: the constructor sets the built-in lines, extendLineRegistry()
    // appends to it in place
    std::vector<LineDef> _registry;

    // registry mapping recombination lines to loaded Case B emissivity tables; filled by
    // initializeRecombinationTables() and extended by extendLineRegistry(), read-only afterwards.
    // Every underlying stab file is bundled into a single dictionary (see StoredTableDictionary)
    // so that opening the full extended catalog does not require a separate memory map -- and
    // thus a separate open file -- per line; StoredTable is move-only, so growing the table vector
    // via push_back()/emplace_back() moves existing elements rather than copying them, which would
    // be unsafe.
    struct RecombRegistry
    {
        bool ready = false;
        StoredTableDictionary<2> dict;
        std::vector<StoredTable<2>> table = std::vector<StoredTable<2>>(numLines);
        std::vector<bool> loaded = std::vector<bool>(numLines, false);
    };
    RecombRegistry _recombRegistry;

    // maps the built-in collisional lines (and any added by extendLineRegistry()) to loaded atomic
    // models and transitions; filled once by initializeAtomicModels(), extended by
    // extendLineRegistry(), read-only otherwise
    struct AtomicLineRegistry
    {
        bool ready = false;
        std::vector<AtomicModel> models;                                   // one per carrier species
        std::vector<std::string> modelNames;                               // species name per model slot
        std::vector<int> lineModel = std::vector<int>(numLines, -1);       // model slot per line, -1 = none
        std::vector<int> lineTransition = std::vector<int>(numLines, -1);  // transition index within the model
    };
    AtomicLineRegistry _atomicRegistry;

    // solves the collisional line luminosity via the atomic-model solver; helper for
    // collisionalLineLuminosity() when an atomic model is available for the line
    double solvedCollisionalLineLuminosity(int lineIdx, double T, double ne, double nIon, double V_cm3) const;
};

//////////////////////////////////////////////////////////////////////

#endif
