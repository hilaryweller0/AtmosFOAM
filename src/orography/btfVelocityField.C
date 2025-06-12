#include "btfVelocityField.H"
#include "btfTransform.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(btfVelocityField, 0);
addToRunTimeSelectionTable(velocityField, btfVelocityField, dict);

btfVelocityField::btfVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
m(mountain::New(dict.subDict("mountain"))),
mode(word(dict.lookup("analyticSolution")) == "horizontalOnly"
        ? analyticSolution::HORIZONTAL_ONLY
        : analyticSolution::VERTICAL_ONLY),
transform(btfTransform(dict))
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
    if (mode == analyticSolution::HORIZONTAL_ONLY)
    {
        return initialHorizontalPositionOf(p, t);
    }
    else
    {
        if (t.value() > 0.0)
        {
            return transform.physicalToComputational(p);
        }
        else
        {
            return p;
        }
    }
}

point btfVelocityField::initialHorizontalPositionOf
(
    const point& p,
    const Time& t
) const
{
    const crossableMountain& mountain = dynamic_cast<const crossableMountain&>(m());

    dimensionedScalar x("x", dimLength, p.x());
    dimensionedScalar t1 = (mountain.start() - x) / u0;
    dimensionedScalar t2 = mountain.timeToCross(u0, H);

    if (t.value() <= t1.value())
    {
        return point((x - t * u0).value(), p.y(), p.z());
    }
    else if (t.value() >= (t1 + t2).value())
    {
        return point((x - (mountain.end() - mountain.start()) - u0 * (t - t2)).value(), p.y(), p.z());
    }
    else
    {
        return p;
    }
}
}
