#include "TestStencils.H"

Foam::localStencil Test::Stencils::twelvePoints()
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
    return localStencil(stencil);
}

Foam::scalarList Test::Coefficients::twelvePoints()
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

Foam::scalarList Test::Coefficients::twoByThreeWithDiagonal()
{
    Foam::scalarList coefficients(6, scalar(0));
    coefficients[0] = 0.5101;
    coefficients[1] = 0.49;
    coefficients[2] = -0.0184;
    coefficients[3] = -0.0155;
    coefficients[4] = 0.0158;
    coefficients[5] = 0.0181;
    return coefficients;
}

Foam::localStencil Test::Stencils::twoPointsInHorizontalLine()
{
    Foam::List<point> stencil(2, point(0, 0, 0));
    stencil[0] = point(-1, 0, 0);
    stencil[1] = point(3, 0, 0);
    return localStencil(stencil);
}

Foam::localStencil Test::Stencils::twoPointsInVerticalLine()
{
    Foam::List<point> stencil(2, point(0, 0, 0));
    stencil[0] = point(0, 1, 0);
    stencil[1] = point(0, -2, 0);
    return localStencil(stencil);
}

Foam::localStencil Test::Stencils::twoByThreeWithDiagonal()
{
    Foam::List<point> stencil(6, point(0, 0, 0));
    stencil[0] = point(-1, 0.0330314, 0);
    stencil[1] = point(0.907926, -2.64e-14, 0);
    stencil[2] = point(-1.13441, 3.05716, 0);
    stencil[3] = point(-0.87524, -2.99221, 0);
    stencil[4] = point(0.907926, 3.02642, 0);
    stencil[5] = point(0.907926, -3.02642, 0);
    return localStencil(stencil);
}

Foam::localStencil Test::Stencils::btfStable()
{
    Foam::List<point> stencil(9, point(0, 0, 0));
    stencil[0] = point(1, -0.0373978, 2.41746e-16);
    stencil[1] = point(-1, 0.0547895, 1.16175e-16);
    stencil[2] = point(2.09757, 2.03321, -5.95858e-16);
    stencil[3] = point(1.9996, -0.0921689, 2.46512e-16);
    stencil[4] = point(1.90163, -2.21755, 6.05113e-16);
    stencil[5] = point(1.05656, 2.08695, -7.20833e-16);
    stencil[6] = point(0.944101, -2.16689, 4.80117e-16);
    stencil[7] = point(-1.02598, 2.18294, 5.03523e-33);
    stencil[8] = point(-0.971564, -2.07859, 5.96625e-16);
    return localStencil(stencil);
}

Foam::localStencil Test::Stencils::slantedCellStable()
{
    Foam::List<point> stencil(14, point(0, 0, 0));
    stencil[0] = point(1, 9.32232e-15, 9.09495e-16);
    stencil[1] = point(-0.292076, 0.555556, 5.25905e-16);
    stencil[2] = point(1, -3.33333, -4.26945e-16);
    stencil[3] = point(1, 3.33333, 8.81693e-16);
    stencil[4] = point(-0.438115, -3.33333, 4.82549e-16);
    stencil[5] = point(-0.438115, 3.7744e-14, 4.83439e-30);
    stencil[6] = point(5, -3.33333, -8.81693e-16);
    stencil[7] = point(5, 1.8872e-14, -4.54747e-16);
    stencil[8] = point(5, 3.33333, -9.37297e-16);
    stencil[9] = point(-0.292076, -3.88889, -4.98103e-16);
    stencil[10] = point(-0.755646, 3.55039, 0);
    stencil[11] = point(3, -3.33333, 9.37297e-16);
    stencil[12] = point(3, 9.32232e-15, 1.19404e-30);
    stencil[13] = point(3, 3.33333, 4.26945e-16);
    return localStencil(stencil);
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
