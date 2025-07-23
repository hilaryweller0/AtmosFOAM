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
radius("radius", dimLength, dict.lookupOrDefault<scalar>("radius", scalar(6.3712e6))),
alpha(dict.lookupOrDefault<scalar>("tilt", scalar(0)))
{};

vector geodesicSolidBodyVelocityField::streamfunctionAt
(
    const point& p,
    scalar time
) const
{
    const scalar u0 = 2 * M_PI * radius.value() / endTime_.value();
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar lon = polarp.lon();

    const scalar psi = - u0 * (Foam::sin(lat) * Foam::cos(alpha) - Foam::cos(lon) * Foam::cos(lat) * Foam::sin(alpha));

    return p/mag(p) * psi * radius.value();
}

point geodesicSolidBodyVelocityField::initialPositionOf
(
    const point& p,
    scalar time
) const
{
    // assume alpha = 0
    
    const scalar u0 = 2 * M_PI * radius.value() / endTime_.value();
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    scalar lon = polarp.lon();

    lon -= u0/radius.value() * time; 

    const polarPoint departureP(lon, lat, radius.value());
    return departureP.cartesian();
}
}
