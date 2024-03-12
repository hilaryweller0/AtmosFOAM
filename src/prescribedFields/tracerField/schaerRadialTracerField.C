#include "schaerRadialTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(schaerRadialTracerField, 0);
addToRunTimeSelectionTable(tracerField, schaerRadialTracerField, dict);

schaerRadialTracerField::schaerRadialTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1.0))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0.0, 9e3))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 1.0, 3e3))),
exponent(dict.lookupOrDefault<scalar>("exponent", scalar(2.0)))
{};

scalar schaerRadialTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar r = mag(vector(d.x()/A.x(), 0, d.z()/A.z()));

    if (r <= 1)
    {
        return rho0*pow(Foam::cos(M_PI*r/2), exponent);
    }
    else
    {
        return 0;
    }
}
}
