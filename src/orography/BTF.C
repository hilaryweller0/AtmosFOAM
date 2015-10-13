#include "BTF.H"
#include "fvCFD.H"

BTF::BTF(const Mountain& mountain, const IOdictionary& dict) :
    mountain(mountain),
    H(readScalar(dict.lookup("domainHeight")))
{};

point BTF::transform(const point& geometric) const
{
    const scalar h = mountain.heightAt(geometric.x());
    const scalar z = geometric.z();
    const scalar z_star = z > H ? z : H * (z - h) / (H - h);
    return point(geometric.x(), geometric.y(), z_star);
}
