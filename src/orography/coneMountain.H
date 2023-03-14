#ifndef CONE_MOUNTAIN_H
#define CONE_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class coneMountain
    : public mountain
{
    public:
    TypeName("coneMountain");
    coneMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;

    private:
    const dimensionedScalar h0_;
    const point centre_;
    const dimensionedScalar radius_;
};

}

#endif

