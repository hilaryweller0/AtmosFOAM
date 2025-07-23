#ifndef SOLID_BODY_ROTATION_VELOCITY_FIELD_H
#define SOLID_BODY_ROTATION_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentVelocityField.H"

namespace Foam
{
class solidBodyRotationVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("solidBodyRotation");
    solidBodyRotationVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, scalar time) const;
    point initialPositionOf(const point& p, scalar time) const;
    const vector rotation_;
    const point centre_;
    const scalar innerRadius_;
    const scalar outerRadius_;
};
}
#endif

