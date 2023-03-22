#include "sphericalLayerTracerField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

defineTypeNameAndDebug(sphericalLayerTracerField, 0);
addToRunTimeSelectionTable(tracerField, sphericalLayerTracerField, dict);

sphericalLayerTracerField::sphericalLayerTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    a_(readScalar(dict.lookup("earthRadius"))),
    z1_(readScalar(dict.lookup("tracerBase"))),
    z2_(readScalar(dict.lookup("tracerTop"))),
    z0_(0.5*(z1_+z2_)),
    H_(z2_ - z1_)
{};

scalar sphericalLayerTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    scalar z = mag(p) - a_;
    scalar q = 0;
    
    if (z > z1_ && z < z2_)
    {
        q = 0.5*(1 + Foam::cos(2*pi*(z-z0_)/H_));
    }

    return q;
}
