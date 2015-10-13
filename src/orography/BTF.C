#include "BTF.H"
#include "fvCFD.H"

BTF::BTF(const Mountain& mountain, const IOdictionary& dict) :
    mountain(mountain),
    H(readScalar(dict.lookup("domainHeight")))
{};

point BTF::transform(const point& geometric) const
{
    return point(geometric.x(), geometric.y(), geometric.z());
}
