#include "radialDistance.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(radialDistance, 0);
addToRunTimeSelectionTable(tracerField, radialDistance, dict);

radialDistance::radialDistance
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3)))
{};

scalar radialDistance::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar r = mag(d);

    return r;
}
