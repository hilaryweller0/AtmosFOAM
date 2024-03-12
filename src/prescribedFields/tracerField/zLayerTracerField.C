#include "zLayerTracerField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
using namespace constant::mathematical;

defineTypeNameAndDebug(zLayerTracerField, 0);
addToRunTimeSelectionTable(tracerField, zLayerTracerField, dict);

zLayerTracerField::zLayerTracerField
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
    z1_(readScalar(dict.lookup("tracerBase"))),
    z2_(readScalar(dict.lookup("tracerTop"))),
    z0_(0.5*(z1_+z2_)),
    H_(z2_ - z1_)
{};

scalar zLayerTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar z = spherical()
             ? mag(p) - earthRadius()
             : p.z();
    scalar q = 0;
    
    if (z > z1_ && z < z2_)
    {
        q = 0.5*(1 + Foam::cos(2*pi*(z-z0_)/H_));
    }

    return q;
}
}
