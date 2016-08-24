#ifndef HORIZONTAL_VELOCITY_FIELD_H
#define HORIZONTAL_VELOCITY_FIELD_H

#include "advectable.H"
#include "velocityField.H"

class horizontalVelocityField final : public velocityField, public advectable
{
    public:
    TypeName("horizontal");
    horizontalVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    const dimensionedScalar u0;
    const dimensionedScalar z1;
    const dimensionedScalar z2;
};

#endif

