#ifndef DOUBLE_COSINE_MOUNTAIN_H
#define DOUBLE_COSINE_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class doubleCosineMountain
    : public mountain
{
    public:
    TypeName("doubleCosineMountain");
    doubleCosineMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;

    private:
    const dimensionedScalar h0_;
    const point centre_;
    const dimensionedScalar radius_;
    const dimensionedScalar h0m_;
    const point centrem_;
    const dimensionedScalar radiusm_;
};

}

#endif

