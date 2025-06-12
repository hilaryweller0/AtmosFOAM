#ifndef GAUSSIAN_MOUNTAIN_H
#define GAUSSIAN_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class gaussianMountain
    : public mountain
{
    public:
    TypeName("gaussianMountain");
    gaussianMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;

    private:
    const dimensionedScalar h0_;
    const point centre_;
    const dimensionedScalar halfwidth_;
};

}

#endif

