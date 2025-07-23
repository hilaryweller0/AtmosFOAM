#ifndef ZERO_VELOCITY_FIELD_H
#define ZERO_VELOCITY_FIELD_H

#include "advectable.H"
#include "divergentVelocityField.H"

// strictly, a zero velocity field is nonDivergent
// but it's easier to specify it without using a streamfunction
namespace Foam
{
class zeroVelocityField final : public divergentVelocityField, public advectable
{
    public:
    TypeName("zeroVelocity");
    zeroVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, scalar time) const;
    point initialPositionOf(const point& p, scalar time) const;
};
}
#endif

