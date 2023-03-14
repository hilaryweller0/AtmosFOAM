#include "cosSqrTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(cosSqrTracerField, 0);
addToRunTimeSelectionTable(tracerField, cosSqrTracerField, dict);

cosSqrTracerField::cosSqrTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    width_(readScalar(dict.lookup("width"))),
    centre_(dict.lookup("centre")),
    maxTracer_(readScalar(dict.lookup("maxTracer"))),
    backgroundTracer_(dict.lookupOrDefault<scalar>("backgroundTracer", 0))
{}

scalar cosSqrTracerField::tracerAt(const point& p, const Time& t) const
{
    if (mag(p - centre_) < width_)
    {
        return backgroundTracer_
             + maxTracer_*0.25*sqr(1+Foam::cos(M_PI*mag(p - centre_)/width_));
    }
    else
    {
        return backgroundTracer_;
    }
}
