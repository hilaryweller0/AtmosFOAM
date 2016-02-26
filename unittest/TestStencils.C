#include "TestStencils.H"

Foam::List<point> twelvePointStencil()
{
    Foam::List<point> stencil(12, point(0, 0, 0));
    stencil[0] = point(-1, 0, 0);
    stencil[1] = point(1, 0, 0);
    stencil[2] = point(-3, 0, 0);
    stencil[3] = point(-5, 0, 0);
    stencil[4] = point(-1, -2, 0);
    stencil[5] = point(1, -2, 0);
    stencil[6] = point(-3, -2, 0);
    stencil[7] = point(-5, -2, 0);
    stencil[8] = point(-1, 2, 0);
    stencil[9] = point(1, 2, 0);
    stencil[10] = point(-3, 2, 0);
    stencil[11] = point(-5, 2, 0);
    return stencil;
}
