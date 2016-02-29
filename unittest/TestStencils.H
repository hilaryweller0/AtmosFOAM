#ifndef TestStencils_H
#define TestStencils_H

#include "fvCFD.H"

Foam::List<point> twelvePointStencil();

Foam::scalarList twelvePointStencilCoefficients();

const scalar twelvePointStencilMatrix[][9] = {
    {1, -1,  0,  1,  0, 0,   -1,   0,   0},
    {1,  1,  0,  1,  0, 0,    1,   0,   0},
    {1, -3,  0,  9,  0, 0,  -27,   0,   0},
    {1, -5,  0, 25,  0, 0, -125,   0,   0},
    {1, -1, -2,  1,  2, 4,   -1,  -2,  -4},
    {1,  1, -2,  1, -2, 4,    1,  -2,   4},
    {1, -3, -2,  9,  6, 4,  -27, -18, -12},
    {1, -5, -2, 25, 10, 4, -125, -50, -20},
    {1, -1,  2,  1, -2, 4,   -1,   2,  -4},
    {1,  1,  2,  1,  2, 4,    1,   2,   4},
    {1, -3,  2,  9, -6, 4,  -27,  18, -12},
    {1, -5,  2, 25,-10, 4, -125,  50, -20}
};

#endif
