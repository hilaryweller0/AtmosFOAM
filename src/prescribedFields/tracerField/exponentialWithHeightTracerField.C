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
    tracerField
    (
        velocityField,
        readBool(dict.lookup("spherical")),
        dict.lookupOrDefault<scalar>("earthRadius", scalar(0))
    ),
    z0_(readScalar(dict.lookup("z0"))),
    T0_(readScalar(dict.lookup("tracerAtz0"))),
    H_(readScalar(dict.lookup("scaleHeight"))),
    offset_(dict.lookupOrDefault<scalar>("offset", scalar(0)))
{};

scalar exponentialWithHeightTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar z = spherical()
             ? mag(p) - earthRadius()
             : p.z();
    z -= z0_;
    return offset_ + T0_*Foam::exp(-z/H_);
}
