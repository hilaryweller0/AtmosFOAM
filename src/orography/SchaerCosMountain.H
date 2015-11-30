#ifndef SCHAER_COS_MOUNTAIN_H
#define SCHAER_COS_MOUNTAIN_H

#include "fvCFD.H"
#include "Mountain.H"

class SchaerCosMountain : public Mountain
{
    public:
    TypeName("SchaerCosMountain");
    SchaerCosMountain(const IOdictionary& dict);
    SchaerCosMountain(const scalar a, const scalar h0, const scalar lambda);
    scalar heightAt(const point& p) const;
    scalar gradientAt(const scalar x) const;
    scalar timeToCross(const scalar u0, const scalar H) const;
    scalar halfWidth() const;
    ~SchaerCosMountain();
    private:
    const scalar a;
    const scalar h0;
    const scalar lambda;
    SchaerSmoothCosMountain *smooth;
    SchaerFineMountain *fine;
};

#endif
