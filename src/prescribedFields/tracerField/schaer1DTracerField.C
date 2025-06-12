#include "schaer1DTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(schaer1DTracerField, 0);
addToRunTimeSelectionTable(tracerField, schaer1DTracerField, dict);

schaer1DTracerField::schaer1DTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 0))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 1, 1)))
{};

scalar schaer1DTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar r = mag(vector(d.x()/A.x(), 0, 0));

    if (r <= 1)
    {
        return rho0*sqr(Foam::cos(M_PI*r/2));
    }
    else
    {
        return 0;
    }
}
}
