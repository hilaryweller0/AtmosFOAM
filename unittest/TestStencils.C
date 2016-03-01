#include "TestStencils.H"

Foam::List<point> Test::Stencils::twelvePoints()
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

Foam::scalarList twelvePointStencilCoefficients()
{
    Foam::scalarList coefficients(12, scalar(0));
    coefficients[0] = 0.875;
    coefficients[1] = 0.297;
    coefficients[2] = -0.141;
    coefficients[3] = -0.031;
    coefficients[4] = 0.031;
    coefficients[5] = 0.008;
    coefficients[6] = -0.086;
    coefficients[7] = 0.047;
    coefficients[8] = 0.031;
    coefficients[9] = 0.008;
    coefficients[10] = -0.086;
    coefficients[11] = 0.047;
    return coefficients;
}

Foam::List<point> Test::Stencils::twoPointsInHorizontalLine()
{
    Foam::List<point> stencil(2, point(0, 0, 0));
    stencil[0] = point(-1, 0, 0);
    stencil[1] = point(3, 0, 0);
    return stencil;
}

Foam::List<point> Test::Stencils::twoPointsInVerticalLine()
{
    Foam::List<point> stencil(2, point(0, 0, 0));
    stencil[0] = point(0, 1, 0);
    stencil[1] = point(0, -2, 0);
    return stencil;
}

const scalar Test::Matrices::twelvePoints[12][9] = {
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

const scalar Test::Matrices::xLinear[2][2] = {
    {1, -1},
    {1, 3}
};

const scalar Test::Matrices::yLinear[2][2] = {
    {1, 1},
    {1, -2}
};
