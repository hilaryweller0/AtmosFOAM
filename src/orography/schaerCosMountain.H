#ifndef SCHAER_COS_MOUNTAIN_H
#define SCHAER_COS_MOUNTAIN_H

#include "mountain.H"

class schaerCosMountain final : public mountain
{
    public:
    TypeName("schaerCos");
    schaerCosMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;

    private:
    const dimensionedScalar h0;
    const scalar a;
    const scalar lambda;
};

#endif

