#ifndef SCHAER_EXP_MOUNTAIN_H
#define SCHAER_EXP_MOUNTAIN_H

#include "dualGradeMountain.H"

namespace Foam
{

class schaerExpMountain final : public dualGradeMountain
{
    public:
    TypeName("schaerExp");
    schaerExpMountain(const dictionary& dict);
    dimensionedScalar coarseHeightAt(const point& p) const;
    dimensionedScalar fineHeightAt(const point& p) const;
    dimensionedScalar start() const;
    dimensionedScalar end() const;

    private:
    const dimensionedScalar h0;
    const scalar a;
    const scalar lambda;
};

}

#endif

