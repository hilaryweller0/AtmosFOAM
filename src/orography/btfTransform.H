#ifndef BTF_TRANSFORM_H
#define BTF_TRANSFORM_H

#include "mountain.H"
#include "terrainFollowingTransform.H"

class btfTransform final : public terrainFollowingTransform
{
    public:
    TypeName("btf");
    btfTransform(const dictionary& dict);
    point physicalToComputational(const point& p) const;
    point computationalToPhysical(const point& p) const;

    private:
    const dimensionedScalar H;
    const autoPtr<mountain> m;
    const scalar groundOffset_;
};

#endif
