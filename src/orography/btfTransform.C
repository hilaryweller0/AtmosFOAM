#include "btfTransform.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(btfTransform, 0);
addToRunTimeSelectionTable(terrainFollowingTransform, btfTransform, dict);

btfTransform::btfTransform(const dictionary& dict)
:
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
m(mountain::New(dict.subDict("mountain")))
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
    const dimensionedScalar h = m->heightAt(p);
    const dimensionedScalar z_star("z_star", dimLength, p.z());
    const dimensionedScalar z = z_star < H ? z_star + h*(1 - z_star/H) : z_star;
    return point(p.x(), p.y(), z.value());
}
