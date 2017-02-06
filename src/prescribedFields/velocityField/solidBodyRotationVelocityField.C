#include "solidBodyRotationVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"

defineTypeNameAndDebug(solidBodyRotationVelocityField, 0);
addToRunTimeSelectionTable(velocityField, solidBodyRotationVelocityField, dict);

solidBodyRotationVelocityField::solidBodyRotationVelocityField(const dictionary& dict)
:
    rotation_(dict.lookup("solidBodyRotation")),
    centre_(dict.lookup("centreOfRotation"))
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

    return magSqr(d)*rotation_;
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
