#include "polynomialTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(polynomialTracerField, 0);
addToRunTimeSelectionTable(tracerField, polynomialTracerField, dict);

polynomialTracerField::polynomialTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    
    a_const(readScalar(dict.lookup("a_const"))),
    
    a_x(dict.lookupOrDefault<scalar>("a_x", scalar(0))),
    a_y(dict.lookupOrDefault<scalar>("a_y", scalar(0))),
    a_z(dict.lookupOrDefault<scalar>("a_z", scalar(0))),
    
    a_xx(dict.lookupOrDefault<scalar>("a_xx", scalar(0))),
    a_xy(dict.lookupOrDefault<scalar>("a_xy", scalar(0))),
    a_xz(dict.lookupOrDefault<scalar>("a_xz", scalar(0))),
    a_yy(dict.lookupOrDefault<scalar>("a_yy", scalar(0))),
    a_yz(dict.lookupOrDefault<scalar>("a_yz", scalar(0))),
    a_zz(dict.lookupOrDefault<scalar>("a_zz", scalar(0))),
    
    a_xxx(dict.lookupOrDefault<scalar>("a_xxx", scalar(0))),
    a_xxy(dict.lookupOrDefault<scalar>("a_xxy", scalar(0))),
    a_xxz(dict.lookupOrDefault<scalar>("a_xxz", scalar(0))),
    a_xyy(dict.lookupOrDefault<scalar>("a_xyy", scalar(0))),
    a_xyz(dict.lookupOrDefault<scalar>("a_xyz", scalar(0))),
    a_xzz(dict.lookupOrDefault<scalar>("a_xzz", scalar(0))),
    a_yyy(dict.lookupOrDefault<scalar>("a_yyy", scalar(0))),
    a_yyz(dict.lookupOrDefault<scalar>("a_yyz", scalar(0))),
    a_yzz(dict.lookupOrDefault<scalar>("a_yzz", scalar(0))),
    a_zzz(dict.lookupOrDefault<scalar>("a_zzz", scalar(0))),
    
    polyMax(dict.lookupOrDefault<scalar>("polyMax", GREAT)),
    polyMin(dict.lookupOrDefault<scalar>("polyMin", -GREAT))
{};

scalar polynomialTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    scalar x = p.x();
    scalar y = p.y();
    scalar z = p.z();
    scalar poly = a_const + a_x*x + a_y*y + a_z*z
         + a_xx*sqr(x) + a_xy*x*y + a_xz*x*z + a_yy*sqr(y) + a_yz*y*z + a_zz*sqr(z)
         + a_xxx*pow(x,3) + a_xxy*sqr(x)*y + a_xxz*sqr(x)*z + a_xyy*x*sqr(y)
         + a_xyz*x*y*z + a_xzz*x*sqr(z)
         + a_yyy*pow(y,3) + a_yyz*sqr(y)*z + a_yzz*y*sqr(z) + a_zzz*pow(z,3);
    return max(min(poly, polyMax), polyMin);
}
}
