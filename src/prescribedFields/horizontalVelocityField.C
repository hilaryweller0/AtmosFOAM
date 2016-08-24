#include "horizontalVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(horizontalVelocityField, 0);
addToRunTimeSelectionTable(velocityField, horizontalVelocityField, dict);

horizontalVelocityField::horizontalVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
z1("zeroVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("zeroVelocityHeight", scalar(4e3))),
z2("maxVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("maxVelocityHeight", scalar(5e3)))
{};

vector horizontalVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    return vector(0,0,0);
}
