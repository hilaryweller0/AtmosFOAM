#ifndef ADVECTABLE_H
#define ADVECTABLE_H

#include "fvCFD.H"

namespace Foam
{

class advectable
{
    public:
    virtual ~advectable() {};
    virtual point initialPositionOf(const point& p, const Time& t) const = 0;
};

}

#endif
