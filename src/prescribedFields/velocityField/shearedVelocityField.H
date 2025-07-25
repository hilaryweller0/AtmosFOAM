#ifndef SHEARED_VELOCITY_FIELD_H
#define SHEARED_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentVelocityField.H"

namespace Foam
{
class shearedVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("sheared");
    shearedVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, scalar time) const;
    point initialPositionOf(const point& p, scalar time) const;
    const vector gradU;
    const vector u0;
    const point p0;
    const scalar acceleration;
    vector gradUhat;
    vector normal;
    const scalar magGradU;
    const scalar magU0;
    vector u0hat;
};
}
#endif

