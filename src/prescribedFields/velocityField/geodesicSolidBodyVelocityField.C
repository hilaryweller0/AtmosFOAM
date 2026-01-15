#include "geodesicSolidBodyVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"

namespace Foam
{
defineTypeNameAndDebug(geodesicSolidBodyVelocityField, 0);
addToRunTimeSelectionTable(velocityField, geodesicSolidBodyVelocityField, dict);

geodesicSolidBodyVelocityField::geodesicSolidBodyVelocityField(const dictionary& dict)
:
nonDivergentVelocityField(dict),
radius
(
    dict.lookupOrDefault<dimensionedScalar>
    (
        "radius", dimensionedScalar(dimLength, scalar(6.3712e6))
    ).value()
),
alpha(dict.lookupOrDefault<scalar>("tilt", scalar(0))),
u0(dict.lookupOrDefault<dimensionedScalar>("u0", 2*M_PI*radius/endTime_).value())
{};

vector geodesicSolidBodyVelocityField::streamfunctionAt
(
    const point& p,
    scalar time
) const
{
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar lon = polarp.lon();

    const scalar psi = - u0 * (Foam::sin(lat) * Foam::cos(alpha) - Foam::cos(lon) * Foam::cos(lat) * Foam::sin(alpha));

    return p/polarp.r() * psi * polarp.r();
}

point geodesicSolidBodyVelocityField::initialPositionOf
(
    const point& p,
    scalar time
) const
{
    // assume alpha = 0
    
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    scalar lon = polarp.lon();

    lon -= u0/polarp.r() * time; 

    const polarPoint departureP(lon, lat, polarp.r());
    return departureP.cartesian();
}
}
