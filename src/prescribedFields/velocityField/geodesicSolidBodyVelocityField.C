#include "geodesicSolidBodyVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"

defineTypeNameAndDebug(geodesicSolidBodyVelocityField, 0);
addToRunTimeSelectionTable(velocityField, geodesicSolidBodyVelocityField, dict);

geodesicSolidBodyVelocityField::geodesicSolidBodyVelocityField(const dictionary& dict)
:
radius("radius", dimLength, dict.lookupOrDefault<scalar>("radius", scalar(6.3712e6))),
alpha(dict.lookupOrDefault<scalar>("tilt", scalar(0))),
endTime("endTime", dimTime, dict.lookupOrDefault<scalar>("endTime", scalar(-1)))
{};

vector geodesicSolidBodyVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    const dimensionedScalar T = (endTime.value() == -1 ) ? t.endTime() : endTime;
    const scalar u0 = 2 * M_PI * radius.value() / T.value();
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar lon = polarp.lon();

    const scalar psi = - u0 * (Foam::sin(lat) * Foam::cos(alpha) - Foam::cos(lon) * Foam::cos(lat) * Foam::sin(alpha));

    return p/mag(p) * psi;
}

point geodesicSolidBodyVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    // assume alpha = 0
    
    const dimensionedScalar T = (endTime.value() == -1 ) ? t.endTime() : endTime;
    const scalar u0 = 2 * M_PI * radius.value() / T.value();
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    scalar lon = polarp.lon();

    lon -= u0/radius.value() * t.value(); 

    const polarPoint departureP(lon, lat, radius.value());
    return departureP.cartesian();
}
