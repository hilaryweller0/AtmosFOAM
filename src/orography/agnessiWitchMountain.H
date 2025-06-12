#ifndef AGNESSIWITCH_MOUNTAIN_H
#define AGNESSIWITCH_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class agnessiWitchMountain
    : public mountain
{
    public:
    TypeName("agnessiWitchMountain");
    agnessiWitchMountain(const dictionary& dict);
    dimensionedScalar heightAt(const point& p) const;

    private:
    const dimensionedScalar h0_;
    const point centre_;
    const dimensionedScalar halfwidth_;
};

}

#endif

