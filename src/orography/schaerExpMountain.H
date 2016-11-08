#ifndef SCHAER_EXP_MOUNTAIN_H
#define SCHAER_EXP_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class schaerExpMountain final : public mountain
{
    public:
    TypeName("schaerExp");
    schaerExpMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;
    dimensionedScalar start() const;
    dimensionedScalar end() const;

    private:
    const dimensionedScalar h0;
    const scalar a;
    const scalar lambda;
};

}

#endif

