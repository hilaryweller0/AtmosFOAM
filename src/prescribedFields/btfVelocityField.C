#include "btfVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(btfVelocityField, 0);
addToRunTimeSelectionTable(velocityField, btfVelocityField, dict);

btfVelocityField::btfVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight")))
{};

vector btfVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    return vector(0, 0, 0);
}

point btfVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    return p;
}
