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

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    const Foam::label faceI = 11;

    Test::PolynomialFit fit(twelvePointStencil(), faceI);

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

TEST_CASE("fit full-size stencil to uniform set of points in local coords")
{
    Test::PolynomialFit fit(twelvePointStencil());

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

// while cubicUpwindCPCFit is a correction on upwind, 
// I've added a test for linear correction to help ensure
// that we don't break those particular code paths
TEST_CASE("fit using linear correction")
{
    bool linearCorrection = true;

    Test::PolynomialFit fit(twelvePointStencil(), linearCorrection);

    check(fit.coefficients(), twelvePointStencilCoefficients());
}

TEST_CASE("dimensions of FixedPolynomialMatrix")
{
    Foam::FixedPolynomialMatrix<cubicUpwindCPCFitPolynomial>
        matrix(twelvePointStencil(), 2);
    scalarRectangularMatrix B = matrix.matrix();     
    CHECK(B.m() == 9);
    CHECK(B.n() == 12);
    CHECK(B[0][0] == approx(1));
    CHECK(B[0][1] == approx(-1));
    CHECK(B[0][2] == approx(0));
    CHECK(B[0][3] == approx(1));
    CHECK(B[0][4] == approx(0));
    CHECK(B[0][5] == approx(0));
    CHECK(B[0][6] == approx(-1));
    CHECK(B[0][7] == approx(0));
    CHECK(B[0][8] == approx(0));
}
