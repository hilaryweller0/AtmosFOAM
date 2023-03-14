#include "schaerRadial3DTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerRadial3DTracerField, 0);
addToRunTimeSelectionTable(tracerField, schaerRadial3DTracerField, dict);

schaerRadial3DTracerField::schaerRadial3DTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 25e3, 3e3)))
{};

scalar schaerRadial3DTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar r = mag(vector(d.x()/A.x(), d.y()/A.y(), d.z()/A.z()));

    if (r <= 1)
    {
        return rho0*sqr(Foam::cos(M_PI*r/2));
    }
    else
    {
        return 0;
    }
}
