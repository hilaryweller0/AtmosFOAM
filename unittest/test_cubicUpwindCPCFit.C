#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "AdaptivePolynomial.H"
#include "cubicUpwindCPCFitPolynomial.H"
#include "TestPolynomialFit.H"
#include "TestStencils.H"
#include "Checks.H"
#include <assert.h>

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    const Foam::label faceI = 11;

    Test::PolynomialFit fit(Test::Stencils::twelvePoints(), faceI);

    check(fit.coefficients(), Test::Coefficients::twelvePoints());
}

TEST_CASE("fit full-size stencil to uniform set of points in local coords")
{
    Test::PolynomialFit fit(Test::Stencils::twelvePoints());

    check(fit.coefficients(), Test::Coefficients::twelvePoints());
}

// while cubicUpwindCPCFit is a correction on upwind, 
// I've added a test for linear correction to help ensure
// that we don't break those particular code paths
TEST_CASE("fit using linear correction")
{
    bool linearCorrection = true;

    Test::PolynomialFit fit(Test::Stencils::twelvePoints(), linearCorrection);

    check(fit.coefficients(), Test::Coefficients::twelvePoints());
}

TEST_CASE("a + bx with two points in horizontal line")
{
    const direction dimensions = 2;
    const localStencil stencil = Test::Stencils::twoPointsInHorizontalLine();
    Foam::AdaptivePolynomial<cubicUpwindCPCFitPolynomial>
        matrix(stencil, dimensions);

    check<2, 2>(matrix.matrix(), Test::Matrices::xLinear);
}

TEST_CASE("a + by with two points in vertical line")
{
    const direction dimensions = 2;
    const localStencil stencil = Test::Stencils::twoPointsInVerticalLine();
    Foam::AdaptivePolynomial<cubicUpwindCPCFitPolynomial>
        matrix(stencil, dimensions);

    check<2, 2>(matrix.matrix(), Test::Matrices::yLinear);
}

/*
 * This test case is taken from a real BTF mesh where the von Neumann
 * analysis detected an instability.  When advectionFoam was run on this mesh,
 * this was the first stencil to become unstable.
 */
TEST_CASE("BTF stable")
{
    Test::PolynomialFit fit(Test::Stencils::btfStable());

    checkStable(fit.coefficients());
}
