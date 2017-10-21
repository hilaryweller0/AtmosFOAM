#include "solidBodyRotationVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"

defineTypeNameAndDebug(solidBodyRotationVelocityField, 0);
addToRunTimeSelectionTable(velocityField, solidBodyRotationVelocityField, dict);

solidBodyRotationVelocityField::solidBodyRotationVelocityField(const dictionary& dict)
:
    rotation_(dict.lookup("solidBodyRotation")),
    centre_(dict.lookup("centreOfRotation")),
    innerRadius_(dict.lookupOrDefault<scalar>("innerRadius", scalar(-1))),
    outerRadius_(dict.lookupOrDefault<scalar>("outerRadius", scalar(-1)))
{};

vector solidBodyRotationVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    // Vector from point p to axis of rotation
    vector d = p - centre_;
    d = d - (d & rotation_)*rotation_/magSqr(rotation_);
    scalar radius = mag(d);
    vector streamFunc = vector::zero;
    
    if (innerRadius_ <=0 || outerRadius_ <= 0 || radius < innerRadius_)
    {
        streamFunc = magSqr(d)*rotation_;
    }
    else if (radius < outerRadius_)
    {
        streamFunc = rotation_*
        (
            sqr(innerRadius_)
          + 2*innerRadius_/(outerRadius_-innerRadius_)*
            (
                outerRadius_*(radius - innerRadius_)
              - 0.5*(sqr(radius) - sqr(innerRadius_))
            )
        );
    }
    else
    {
        streamFunc = rotation_*
        (
            sqr(innerRadius_)
          + 2*innerRadius_/(outerRadius_-innerRadius_)*
            (
                outerRadius_*(outerRadius_-innerRadius_)
              - 0.5*(sqr(outerRadius_) - sqr(innerRadius_))
            )
        );
    }

    return streamFunc;
}

point solidBodyRotationVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    // rotation matrix
    tensor R = Ra(rotation_, -mag(rotation_)*t.value());
    
    return (R & (p - centre_)) + centre_;
}
