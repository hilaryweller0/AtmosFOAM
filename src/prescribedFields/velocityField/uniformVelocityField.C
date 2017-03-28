#include "uniformVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"
#include "sphericalVector.H"

defineTypeNameAndDebug(uniformVelocityField, 0);
addToRunTimeSelectionTable(velocityField, uniformVelocityField, dict);

uniformVelocityField::uniformVelocityField(const dictionary& dict)
:
v(dict.lookupOrDefault<vector>("velocity", vector(1,1,1)))
{};

vector uniformVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    return v;
}

point uniformVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    scalar time = t.value();
    return point
    (
        p.x() - v[0] * time,
        p.y() - v[1] * time,
        p.z() - v[2] * time
    );
}
