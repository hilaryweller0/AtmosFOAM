#ifndef TestStencils_H
#define TestStencils_H

#include "fvCFD.H"
#include "localStencil.H"

// TODO: I wonder whether each test should have an object with a stencil, expected matrix, expected coefficients etc
namespace Test {
    class Stencils
    {
    public:
        static Foam::localStencil twelvePoints();
        static Foam::localStencil twoPointsInHorizontalLine();
        static Foam::localStencil twoPointsInVerticalLine();
        static Foam::localStencil twoByThreeWithDiagonal();
        static Foam::localStencil btfStable();
        static Foam::localStencil slantedCellStable();
    };

    class Coefficients
    {
    public:
        static Foam::scalarList twelvePoints();
        static Foam::scalarList twoByThreeWithDiagonal();
    };

    class Matrices
    {
    public:
        static const scalar twelvePoints[12][9];
        static const scalar xLinear[2][2];
        static const scalar yLinear[2][2];
    };
}

#endif
