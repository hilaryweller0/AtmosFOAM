#include "btfTransform.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(btfTransform, 0);
addToRunTimeSelectionTable(terrainFollowingTransform, btfTransform, dict);

btfTransform::btfTransform(const dictionary& dict)
:
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
m(mountain::New(dict.subDict("mountain"))),
groundOffset_(dict.lookupOrDefault<scalar>("groundOffset", scalar(0)))
{}

point btfTransform::physicalToComputational(const point& p) const
{
    const dimensionedScalar h = m->heightAt(p);
    const scalar z = p.z();
    const scalar z_star = z > H.value() ? z : (H * (z - h.value()) / (H - h)).value();
    return point(p.x(), p.y(), z_star);
}

point btfTransform::computationalToPhysical(const point& p) const
{
    scalar h = m->heightAt(p).value();
    scalar z_star = p.z();
    scalar z = z_star < H.value() ? z_star + h*(1 - z_star/H.value())
                                        : z_star;
    // Apply offset
    if (z_star < mag(groundOffset_))
    {
        z += groundOffset_;
    }
    return point(p.x(), p.y(), z);
}
