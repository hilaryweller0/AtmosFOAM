#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "FixedPolynomialMatrix.H"
#include "cubicUpwindCPCFitPolynomial.H"
#include "TestPolynomialFit.H"
#include "TestStencils.H"
#include "Checks.H"
#include <assert.h>

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    const Foam::label faceI = 11;

    Test::PolynomialFit fit(Test::Stencils::twelvePoints(), faceI);

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

TEST_CASE("fit full-size stencil to uniform set of points in local coords")
{
    Test::PolynomialFit fit(Test::Stencils::twelvePoints());

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

// while cubicUpwindCPCFit is a correction on upwind, 
// I've added a test for linear correction to help ensure
// that we don't break those particular code paths
TEST_CASE("fit using linear correction")
{
    bool linearCorrection = true;

    Test::PolynomialFit fit(Test::Stencils::twelvePoints(), linearCorrection);

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

TEST_CASE("12x9 FixedPolynomialMatrix")
{
    const direction dimensions = 2;
    Foam::FixedPolynomialMatrix<cubicUpwindCPCFitPolynomial>
        matrix(Test::Stencils::twelvePoints(), dimensions);
    scalarRectangularMatrix B = matrix.matrix();

    const scalar (&expected)[12][9] = twelvePointStencilMatrix;

    CHECK(B.n() == 12);
    CHECK(B.m() == 9);

    for (int i=0; i<B.n(); i++)
    {
        for (int j=0; j<B.m(); j++)
        {
            CHECK(B[i][j] == approx(expected[i][j]));
        }
    }
}

TEST_CASE("a + bx with two points in horizontal line")
{
    const direction dimensions = 2;
    Foam::FixedPolynomialMatrix<cubicUpwindCPCFitPolynomial>
        matrix(Test::Stencils::twoPointsInHorizontalLine(), dimensions);
}
