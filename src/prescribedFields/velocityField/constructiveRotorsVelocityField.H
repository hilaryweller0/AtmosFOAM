#ifndef CONSTRUVTIVE_ROTORS_VELOCITY_FIELD_H
#define CONSTRUCTIVE_ROTORS_VELOCITY_FIELD_H

#include "advectable.H"
#include "divergentVelocityField.H"

// strictly, a uniform velocity field is nonDivergent
// but it's easier to specify it without using a streamfunction
class constructiveRotorsVelocityField final : public divergentVelocityField, public advectable
{
    public:
    TypeName("constructiveRotors");
    constructiveRotorsVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    const scalar vmax;
    const scalar xmax;
    const scalar ymax;
    const point center;
};

#endif


