/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "gasContinuumEmissionSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void gasContinuumEmissionSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, filename(), "lambda(m),logU(1),Z(1)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> gasContinuumEmissionSEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>{
        {"logU"}, {"metallicity"},{"IonisingLum", "Luminosity", "W"},{"EmissionBool"}
    };
}

////////////////////////////////////////////////////////////////////

// Incomplete: Throw error if outside range

Range gasContinuumEmissionSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}


////////////////////////////////////////////////////////////////////


double gasContinuumEmissionSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{

    double logU         = parameters[0];
    double Z            = parameters[1];
    double IonisingLum  = parameters[2];
    double EmissionBool = parameters[3];

    return IonisingLum * EmissionBool * _table(wavelength, logU, Z);
}

////////////////////////////////////////////////////////////////////

double gasContinuumEmissionSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                    const Array& parameters) const
{
    double logU         = parameters[0];
    double Z            = parameters[1];
    double IonisingLum  = parameters[2];
    double EmissionBool = parameters[3];

    return IonisingLum * EmissionBool * _table.cdf(lambdav, pv, Pv, wavelengthRange, logU, Z);
}

////////////////////////////////////////////////////////////////////

 
