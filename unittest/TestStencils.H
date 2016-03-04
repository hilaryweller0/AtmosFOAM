#ifndef TestStencils_H
#define TestStencils_H

#include "fvCFD.H"

// TODO: I wonder whether each test should have an object with a stencil, expected matrix, expected coefficients etc
namespace Test {
    class Stencils
    {
    public:
        static Foam::List<point> twelvePoints();
        static Foam::List<point> twoPointsInHorizontalLine();
        static Foam::List<point> twoPointsInVerticalLine();
    };

    class Matrices
    {
    public:
        static const scalar twelvePoints[12][9];
        static const scalar xLinear[2][2];
        static const scalar yLinear[2][2];
    };
}


Foam::scalarList twelvePointStencilCoefficients();

#endif
