#include "geodesicSinusoidalTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(geodesicSinusoidalTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicSinusoidalTracerField, dict);

geodesicSinusoidalTracerField::geodesicSinusoidalTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    Rsphere(readScalar(dict.lookup("Rsphere"))),
    c(dict.lookupOrDefault<scalar>("c", scalar(0))),
    
    a_x(dict.lookupOrDefault<scalar>("a_x", scalar(1))),
    a_y(dict.lookupOrDefault<scalar>("a_y", scalar(1))),
    a_z(dict.lookupOrDefault<scalar>("a_z", scalar(1))),
    
    b_x(dict.lookupOrDefault<scalar>("b_x", scalar(0))),
    b_y(dict.lookupOrDefault<scalar>("b_y", scalar(0))),
    b_z(dict.lookupOrDefault<scalar>("b_z", scalar(0))),
    
    k_x(dict.lookupOrDefault<scalar>("k_x", scalar(1))),
    k_y(dict.lookupOrDefault<scalar>("k_y", scalar(1))),
    k_z(dict.lookupOrDefault<scalar>("k_z", scalar(1))),
    
    L_z(dict.lookupOrDefault<scalar>("L_z", scalar(1))),
    
    phi_x(dict.lookupOrDefault<scalar>("phi_x", scalar(0))),
    phi_y(dict.lookupOrDefault<scalar>("phi_y", scalar(0))),
    phi_z(dict.lookupOrDefault<scalar>("phi_z", scalar(0)))
{};

scalar geodesicSinusoidalTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar x = Foam::atan2(p.y(), p.x());
    scalar y = Foam::asin(p.z()/mag(p));
    scalar z = mag(p) - Rsphere;

    return c + (a_x + b_x*Foam::cos(k_x*(x-phi_x)))
             * (a_y + b_y*Foam::cos(k_y*(y-phi_y)))
             * (a_z + b_z*Foam::cos(2*M_PI*k_z*(z-phi_z)/L_z));
}
}
