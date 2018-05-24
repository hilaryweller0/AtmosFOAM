#include "cosineTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(cosineTracerField, 0);
addToRunTimeSelectionTable(tracerField, cosineTracerField, dict);

cosineTracerField::cosineTracerField
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

scalar cosineTracerField::tracerAt(const point& p, const Time& t) const
{
    if (mag(p - centre_) < (width_*2))
    {
        return maxTracer_*0.5*(1+Foam::cos(M_PI*mag(p - centre_)/(width_*2)));
    }
    else
    {
        return 0;
    }
}
