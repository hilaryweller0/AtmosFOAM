#ifndef GEODESIC_SOLID_BODY_VELOCITY_FIELD_H
#define GEODESIC_SOLID_BODY_VELOCITY_FIELD_H

#include "velocityField.H"

class geodesicSolidBodyVelocityField final : public velocityField
{
    public:
    TypeName("geodesicSolidBody");
    geodesicSolidBodyVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    const dimensionedScalar radius;
    const scalar alpha;
};

#endif

