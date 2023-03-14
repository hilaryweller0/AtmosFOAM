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
    width_(readScalar(dict.lookup("width"))),
    centre_(dict.lookup("centre")),
    maxTracer_(readScalar(dict.lookup("maxTracer")))
{}

scalar coneTracerField::tracerAt(const point& p, const Time& t) const
{
    return max(maxTracer_*(1-mag(p - centre_)/width_), 0);
}

