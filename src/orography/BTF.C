#include "BTF.H"
#include "fvCFD.H"

BTF::BTF(const autoPtr<Mountain> mountain, const IOdictionary& dict) :
    mountain(mountain),
    H(readScalar(dict.lookup("domainHeight")))
{};

point BTF::computationalToPhysical(const point& computational) const
{
    const scalar h = mountain->heightAt(computational);
    const scalar z_star = computational.z();
    const scalar z = (z_star < H) ? z_star + h*(1 - z_star/H) : z_star;
    return point(computational.x(), computational.y(), z);
}

point BTF::physicalToComputational(const point& geometric) const
{
    const scalar h = mountain->heightAt(geometric);
    const scalar z = geometric.z();
    const scalar z_star = z > H ? z : H * (z - h) / (H - h);
    return point(geometric.x(), geometric.y(), z_star);
}
