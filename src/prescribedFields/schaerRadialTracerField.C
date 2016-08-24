#include "schaerRadialTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerRadialTracerField, 0);
addToRunTimeSelectionTable(tracerField, schaerRadialTracerField, dict);

schaerRadialTracerField::schaerRadialTracerField(const dictionary& dict)
:
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 0, 3e3)))
{};

scalar schaerRadialTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    return 0;
}
