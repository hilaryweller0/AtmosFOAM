#include "zeroVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(zeroVelocityField, 0);
addToRunTimeSelectionTable(velocityField, zeroVelocityField, dict);

zeroVelocityField::zeroVelocityField(const dictionary& dict)
{};

vector zeroVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    return vector::zero;
}

point zeroVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    return p;
}

}
