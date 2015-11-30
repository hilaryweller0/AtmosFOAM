#ifndef SCHAER_EXP_MOUNTAIN_H
#define SCHAER_EXP_MOUNTAIN_H

#include "fvCFD.H"
#include "Mountain.H"

class SchaerExpMountain : public Mountain
{
    public:
    TypeName("SchaerExpMountain");
    SchaerExpMountain(const IOdictionary& dict);
    SchaerExpMountain(const scalar a, const scalar h0, const scalar lambda);
    scalar heightAt(const point& p) const;
    scalar gradientAt(const scalar x) const;
    scalar timeToCross(const scalar u0, const scalar H) const;
    scalar halfWidth() const;
    ~SchaerExpMountain();
    private:
    const scalar a;
    const scalar h0;
    const scalar lambda;
    SchaerSmoothExpMountain *smooth;
    SchaerFineMountain *fine;
};

#endif
