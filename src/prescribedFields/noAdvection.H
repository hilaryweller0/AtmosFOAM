#ifndef NO_ADVECTION_H
#define NO_ADVECTION_H

#include "advectable.H"

namespace Foam
{

class noAdvection : public advectable
{
    public:
        point initialPositionOf(const point& p, scalar time) const;
        point initialPositionOf(const point& p, const Time& t) const
            { return initialPositionOf(p, t.value()); }
};

}

#endif

