#ifndef ADVECTABLE_H
#define ADVECTABLE_H

#include "fvMesh.H"

namespace Foam
{

class advectable
{
    public:
    virtual ~advectable() {};
    virtual point initialPositionOf(const point& p, scalar time) const = 0;
    point initialPositionOf(const point& p, const Time& t) const
        { return initialPositionOf(p, t.value()); }
};

}

#endif
