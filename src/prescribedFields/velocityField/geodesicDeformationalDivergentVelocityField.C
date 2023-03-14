#include "geodesicDeformationalDivergentVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"
#include "sphericalVector.H"

defineTypeNameAndDebug(geodesicDeformationalDivergentVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    geodesicDeformationalDivergentVelocityField,
    dict
);

geodesicDeformationalDivergentVelocityField::geodesicDeformationalDivergentVelocityField
(
    const dictionary& dict
)
:
    radius("radius", dimLength, readScalar(dict.lookup("radius")))
{};

vector geodesicDeformationalDivergentVelocityField::velocityAt
(
        const point& p,
        const Time& t
) const
{
    const dimensionedScalar T = t.endTime();

    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar lon = polarp.lon();

    const dimensionedScalar lonPrime = lon - 2 * M_PI * t / T;
    
    dimensionedScalar u = -5 * radius / T * sqr(Foam::sin(lonPrime/2)) * Foam::sin(2*lat) * sqr(Foam::cos(lat)) * Foam::cos(M_PI * t / T) + 2*M_PI*radius/T * Foam::cos(lat);
    dimensionedScalar v = 5/2 * radius / T * Foam::sin(lonPrime) * pow3(Foam::cos(lat)) * Foam::cos(M_PI * t / T);

    sphericalVector localWind(u.value(), v.value(), 0);
    sphericalVector sphericalp(p);

    return localWind.toCartesian(sphericalp);
}

