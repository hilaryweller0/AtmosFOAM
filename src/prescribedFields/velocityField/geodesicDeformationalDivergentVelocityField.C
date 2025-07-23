#include "geodesicDeformationalDivergentVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"
#include "sphericalVector.H"

namespace Foam
{
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
    geodesicVelocityField(dict)
{};

vector geodesicDeformationalDivergentVelocityField::velocityAt
(
    const point& p,
    scalar time
) const
{
    scalar a = earthRadius_.value();
    const scalar T = endTime_.value();

    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar lon = polarp.lon();

    const scalar lonPrime = lon - 2 * M_PI * time / T;
    
    dimensionedScalar u = -5*a/endTime_*sqr(Foam::sin(lonPrime/2))
            *Foam::sin(2*lat)*sqr(Foam::cos(lat))*Foam::cos(M_PI*time/T)
          + 2*M_PI*a/endTime_*Foam::cos(lat);
    dimensionedScalar v = 5/2*a/endTime_*Foam::sin(lonPrime)*pow3(Foam::cos(lat))
         *Foam::cos(M_PI*time/T);

    sphericalVector localWind(u.value(), v.value(), 0);
    sphericalVector sphericalp(p);

    return localWind.toCartesian(sphericalp);
}
}
