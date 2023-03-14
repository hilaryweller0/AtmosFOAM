#ifndef CROSSABLE_MOUNTAIN_H
#define CROSSABLE_MOUNTAIN_H

#include "fvCFD.H"
#include "mountain.H"

namespace Foam
{

class crossableMountain : public mountain
{
    public:
    virtual dimensionedScalar start() const = 0;
    virtual dimensionedScalar end() const = 0;
    virtual dimensionedScalar timeToCross
    (
        const dimensionedScalar u0,
        const dimensionedScalar H
    ) const = 0;
};

}

#endif


