#include "uniformVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(uniformVelocityField, 0);
addToRunTimeSelectionTable(velocityField, uniformVelocityField, dict);

uniformVelocityField::uniformVelocityField(const dictionary& dict)
:
    divergentVelocityField(dict),
    v(dict.lookup("velocity")),
    acceleration(dict.lookupOrDefault<vector>("acceleration", vector::zero))
{};

vector uniformVelocityField::velocityAt
(
    const point& p,
    scalar time
) const
{
    return v + acceleration*time;
}

point uniformVelocityField::initialPositionOf
(
    const point& p,
    scalar time
) const
{
    vector dist = v*time + 0.5*acceleration*sqr(time);
    return point(p - dist);
}

}
