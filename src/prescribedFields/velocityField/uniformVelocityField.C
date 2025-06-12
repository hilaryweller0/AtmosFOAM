#include "uniformVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(uniformVelocityField, 0);
addToRunTimeSelectionTable(velocityField, uniformVelocityField, dict);

uniformVelocityField::uniformVelocityField(const dictionary& dict)
:
    v(dict.lookup("velocity")),
    acceleration(dict.lookupOrDefault<vector>("acceleration", vector::zero))
{};

vector uniformVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    return v + acceleration*t.value();
}

point uniformVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    scalar time = t.value();
    vector dist = v*time + 0.5*acceleration*sqr(time);
    return point(p - dist);
}

}
