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
    const scalar a;
    const scalar h0;
    const scalar lambda;
};

#endif

