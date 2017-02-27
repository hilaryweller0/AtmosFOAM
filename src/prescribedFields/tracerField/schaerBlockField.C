#include "schaerBlockField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerBlockField, 0);
addToRunTimeSelectionTable(tracerField, schaerBlockField, dict);

schaerBlockField::schaerBlockField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 1, 3e3)))
{};

scalar schaerBlockField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar r = mag(vector(d.x()/A.x(), 0, d.z()/A.z()));

    if (r <= 1)
    {
        return rho0;
    }
    else
    {
        return 0;
    }
}
