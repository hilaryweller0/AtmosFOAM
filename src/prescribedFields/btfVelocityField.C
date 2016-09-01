#include "btfVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(btfVelocityField, 0);
addToRunTimeSelectionTable(velocityField, btfVelocityField, dict);

btfVelocityField::btfVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
m(autoPtr<crossableMountain>(dynamic_cast<crossableMountain*>(mountain::New(dict.subDict("mountain")).ptr())))
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
    dimensionedScalar x("x", dimLength, p.x());
    dimensionedScalar t1 = (m->start() - x) / u0;
    dimensionedScalar t2 = m->timeToCross(u0, H);

    if (t.value() <= t1.value())
    {
        return point((x - t * u0).value(), p.y(), p.z());
    }
    else if (t.value() >= (t1 + t2).value())
    {
        dimensionedScalar t3 = t - t2 - t1;
        return point((x - (t1 + t3) * u0 - (m->end() - m->start())).value(), p.y(), p.z());
    }
    else
    {
        return p;
    }
}
