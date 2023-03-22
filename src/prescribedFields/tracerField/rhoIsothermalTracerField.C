#include "rhoIsothermalTracerField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

defineTypeNameAndDebug(rhoIsothermalTracerField, 0);
addToRunTimeSelectionTable(tracerField, rhoIsothermalTracerField, dict);

rhoIsothermalTracerField::rhoIsothermalTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    a_(readScalar(dict.lookup("earthRadius"))),
    rho0_(readScalar(dict.lookup("rho0"))),
    H_(readScalar(dict.lookup("scaleHeight")))
{};

scalar rhoIsothermalTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar z = mag(p) - a_;
    return rho0_*Foam::exp(-z/H_);
}
