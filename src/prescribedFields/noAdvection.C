#include "noAdvection.H"

point noAdvection::initialPositionOf
(
        const point& p,
        const Time& t
) const
{
    return p;
}
