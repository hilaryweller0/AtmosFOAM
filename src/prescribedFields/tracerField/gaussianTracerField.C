#include "gaussianTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(gaussianTracerField, 0);
addToRunTimeSelectionTable(tracerField, gaussianTracerField, dict);

gaussianTracerField::gaussianTracerField
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

scalar gaussianTracerField::tracerAt(const point& p, const Time& t) const
{
    return maxTracer_*Foam::exp(-0.5*(magSqr(p - centre_)/sqr(width_)));
}

}
