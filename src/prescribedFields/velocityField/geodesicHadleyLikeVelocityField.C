#include "geodesicHadleyLikeVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

#include "polarPoint.H"
#include "sphericalVector.H"

namespace Foam
{
using namespace constant::mathematical;

defineTypeNameAndDebug(geodesicHadleyLikeVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    geodesicHadleyLikeVelocityField,
    dict
);

geodesicHadleyLikeVelocityField::
geodesicHadleyLikeVelocityField(const dictionary& dict)
:
    geodesicVelocityField(dict),
    H_(readScalar(dict.lookup("scaleHeight"))),
    ztop_(readScalar(dict.lookup("ztop"))),
    u0_(readScalar(dict.lookup("u0"))),
    w0_(readScalar(dict.lookup("w0"))),
    tau_(readScalar(dict.lookup("endTime"))),
    rho0_(readScalar(dict.lookup("rho0"))),
    K_(readLabel(dict.lookup("nOverturningCells")))
{};

vector geodesicHadleyLikeVelocityField::velocityAt
(
    const point& p,
    scalar time
) const
{
    scalar a = geodesicVelocityField::earthRadius_.value();
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar z = polarp.r() - a;

    const scalar rho = rho0_*Foam::exp(-z/H_);

    sphericalVector localWind
    (
        rho*u0_*Foam::cos(lat),
        -a*w0_*pi/(K_*ztop_)*Foam::cos(lat)*Foam::sin(K_*lat)
            *Foam::cos(pi*z/ztop_)*Foam::cos(pi*time/tau_),
        w0_/K_*
        (
            -2*Foam::sin(K_*lat)*Foam::sin(lat)
           + K_*Foam::cos(lat)*Foam::cos(K_*lat)
        )*Foam::sin(pi*z/ztop_)*Foam::cos(pi*time/tau_)
    );

    sphericalVector sphericalp(p);
    return localWind.toCartesian(sphericalp);
}
}
