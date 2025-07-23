#include "noAdvection.H"

namespace Foam
{

point noAdvection::initialPositionOf(const point& p, scalar time) const
{
    return p;
}

}
