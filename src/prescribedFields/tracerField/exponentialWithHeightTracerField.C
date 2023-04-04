#include "exponentialWithHeightTracerField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

defineTypeNameAndDebug(exponentialWithHeightTracerField, 0);
addToRunTimeSelectionTable(tracerField, exponentialWithHeightTracerField, dict);

exponentialWithHeightTracerField::exponentialWithHeightTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    z0_(readScalar(dict.lookup("z0"))),
    T0_(readScalar(dict.lookup("tracerAtz0"))),
    H_(readScalar(dict.lookup("scaleHeight")))
{};

scalar exponentialWithHeightTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar z = p.z() - z0_;
    return T0_*Foam::exp(-z/H_);
}
