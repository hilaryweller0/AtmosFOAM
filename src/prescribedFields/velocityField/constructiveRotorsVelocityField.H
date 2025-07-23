#ifndef CONSTRUVTIVE_ROTORS_VELOCITY_FIELD_H
#define CONSTRUCTIVE_ROTORS_VELOCITY_FIELD_H

#include "advectable.H"
#include "nonDivergentVelocityField.H"

// strictly, a uniform velocity field is nonDivergent
// but it's easier to specify it without using a streamfunction
namespace Foam
{
class constructiveRotorsVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("constructiveRotors");
    constructiveRotorsVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, scalar time) const;
    point initialPositionOf(const point& p, scalar time) const;
    const scalar vmax;
    const scalar xmax;
    const scalar ymax;
    const point center;
};
}
#endif

