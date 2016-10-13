#ifndef GEODESIC_SOLID_BODY_VELOCITY_FIELD_H
#define GEODESIC_SOLID_BODY_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentSphericalVelocityField.H"

class geodesicSolidBodyVelocityField final : public nonDivergentSphericalVelocityField, public advectable
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

#endif

