#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "AdaptivePolynomial.H"
#include "FixedPolynomial.H"
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

TEST_CASE("12x9 FixedPolynomial")
{
    const direction dimensions = 2;
    const List<point> stencil = Test::Stencils::twelvePoints();
    Foam::FixedPolynomial<cubicUpwindCPCFitPolynomial>
        matrix(stencil, dimensions);

    check<12, 9>(matrix.matrix(), Test::Matrices::twelvePoints);
}

TEST_CASE("a + bx with two points in horizontal line")
{
    const direction dimensions = 2;
    const List<point> stencil = Test::Stencils::twoPointsInHorizontalLine();
    Foam::AdaptivePolynomial<cubicUpwindCPCFitPolynomial>
        matrix(stencil, dimensions);

    check<2, 2>(matrix.matrix(), Test::Matrices::xLinear);
}

TEST_CASE("a + by with two points in vertical line")
{
    const direction dimensions = 2;
    const List<point> stencil = Test::Stencils::twoPointsInVerticalLine();
    Foam::AdaptivePolynomial<cubicUpwindCPCFitPolynomial>
        matrix(stencil, dimensions);

    check<2, 2>(matrix.matrix(), Test::Matrices::yLinear);
}
