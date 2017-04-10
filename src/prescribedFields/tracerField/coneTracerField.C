#include "coneTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(coneTracerField, 0);
addToRunTimeSelectionTable(tracerField, coneTracerField, dict);

coneTracerField::coneTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    width_(dict.lookupOrDefault<scalar>("width", scalar(8))),
    centre_(dict.lookupOrDefault<point>("centre", point(100, 0, 100))),
    maxTracer_(dict.lookupOrDefault<scalar>("maxTracer", scalar(1)))
{}

scalar coneTracerField::tracerAt(const point& p, const Time& t) const
{
    return max(maxTracer_*(1-mag(p - centre_)/width_), 0);
}

