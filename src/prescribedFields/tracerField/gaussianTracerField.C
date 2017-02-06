#include "gaussianTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(gaussianTracerField, 0);
addToRunTimeSelectionTable(tracerField, gaussianTracerField, dict);

gaussianTracerField::gaussianTracerField
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

scalar gaussianTracerField::tracerAt(const point& p, const Time& t) const
{
    return maxTracer_*Foam::exp(-0.5*(magSqr(p - centre_)/sqr(width_)));
}

