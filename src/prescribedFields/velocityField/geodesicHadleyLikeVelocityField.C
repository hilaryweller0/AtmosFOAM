#include "geodesicHadleyLikeVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

#include "polarPoint.H"
#include "sphericalVector.H"

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
    a_(readScalar(dict.lookup("earthRadius"))),
    H_(readScalar(dict.lookup("scaleHeight"))),
    ztop_(readScalar(dict.lookup("ztop"))),
    u0_(readScalar(dict.lookup("u0"))),
    w0_(readScalar(dict.lookup("w0"))),
    tau_(readScalar(dict.lookup("endTime"))),
    K_(readLabel(dict.lookup("nOverturningCells")))
{};

vector geodesicHadleyLikeVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    const polarPoint& polarp = convertToPolar(p);
    const scalar lat = polarp.lat();
    const scalar z = polarp.r() - a_;

    const scalar rho = Foam::exp(-z/H_);

    sphericalVector localWind
    (
        u0_*Foam::cos(lat),
        -a_*w0_*pi/(K_*ztop_*rho)*Foam::cos(lat)*Foam::sin(K_*lat)
            *Foam::cos(pi*z/ztop_)*Foam::cos(pi*t.value()/tau_),
        w0_/(K_*rho)*
        (
            -2*Foam::sin(K_*lat)*Foam::sin(lat)
           + K_*Foam::cos(lat)*Foam::cos(K_*lat)
        )*Foam::sin(pi*z/ztop_)*Foam::cos(pi*t.value()/tau_)
    );

    sphericalVector sphericalp(p);
    return localWind.toCartesian(sphericalp);
}
