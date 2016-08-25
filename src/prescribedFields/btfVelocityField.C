#include "btfVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(btfVelocityField, 0);
addToRunTimeSelectionTable(velocityField, btfVelocityField, dict);

btfVelocityField::btfVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
m(mountain::New(dict.subDict("mountain")))
{};

vector btfVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    const vector unitNormal(0, -1, 0);
    if (p.z() > H.value())
    {
        return unitNormal * -u0.value() * p.z();
    }
    else
    {
        dimensionedScalar h = m->heightAt(p);
        return unitNormal * (-u0 * (H * (p.z() - h.value()) / (H - h))).value();
    }
}

point btfVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    return p;
}
