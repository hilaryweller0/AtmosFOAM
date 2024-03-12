#include "cosSqrTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(cosSqrTracerField, 0);
addToRunTimeSelectionTable(tracerField, cosSqrTracerField, dict);

cosSqrTracerField::cosSqrTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    halfWidth_(readScalar(dict.lookup("halfWidth"))),
    centre_(dict.lookup("centre")),
    maxTracer_(readScalar(dict.lookup("maxTracer"))),
    backgroundTracer_(readScalar(dict.lookup("backgroundTracer")))
{}

scalar cosSqrTracerField::tracerAt(const point& p, const Time& t) const
{
    if (mag(p - centre_) < halfWidth_)
    {
        return backgroundTracer_
             + maxTracer_*sqr(Foam::cos(0.5*M_PI*mag(p - centre_)/halfWidth_));
    }
    else
    {
        return backgroundTracer_;
    }
}

}
