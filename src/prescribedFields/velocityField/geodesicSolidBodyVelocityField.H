#ifndef GEODESIC_SOLID_BODY_VELOCITY_FIELD_H
#define GEODESIC_SOLID_BODY_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentVelocityField.H"

namespace Foam
{
class geodesicSolidBodyVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("geodesicSolidBody");
    geodesicSolidBodyVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    const dimensionedScalar radius;
    const scalar alpha;
    dimensionedScalar endTime;
};
}
#endif

