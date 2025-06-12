#include "noAdvection.H"

namespace Foam
{

point noAdvection::initialPositionOf
(
        const point& p,
        const Time& t
) const
{
    return p;
}

}
