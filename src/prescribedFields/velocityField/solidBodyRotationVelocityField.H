#ifndef SOLID_BODY_ROTATION_VELOCITY_FIELD_H
#define SOLID_BODY_ROTATION_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentVelocityField.H"

class solidBodyRotationVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("solidBodyRotation");
    solidBodyRotationVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    const dimensionedScalar u0;
    const dimensionedScalar z1;
    const dimensionedScalar z2;
};

#endif

