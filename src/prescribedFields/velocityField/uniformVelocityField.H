#ifndef UNIFORM_VELOCITY_FIELD_H
#define UNIFORM_VELOCITY_FIELD_H

#include "advectable.H"
#include "divergentVelocityField.H"

// strictly, a uniform velocity field is nonDivergent
// but it's easier to specify it without using a streamfunction
namespace Foam
{
class uniformVelocityField final : public divergentVelocityField, public advectable
{
    public:
    TypeName("uniform");
    uniformVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    const vector v;
    const vector acceleration;
};
}
#endif

