#ifndef DOUBLE_CYLINDER_MOUNTAIN_H
#define DOUBLE_CYLINDER_MOUNTAIN_H

#include "mountain.H"

namespace Foam
{

class doubleCylinderMountain
    : public mountain
{
    public:
    TypeName("doubleCylinderMountain");
    doubleCylinderMountain(const dictionary& dict);
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

