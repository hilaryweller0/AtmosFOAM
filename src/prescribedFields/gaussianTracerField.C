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
dxMid("width", dimLength, dict.lookupOrDefault<scalar>("width", scalar(8))),
centre(dict.lookupOrDefault<point>("centre", point(100, 0, 100)))
{};

scalar gaussianTracerField::tracerAt(const point& p, const Time& t) const
{
    const point diff = p - centre;
    return Foam::exp(- sqr(diff.x()/dxMid.value()) - sqr(diff.z()/dxMid.value()));
}

