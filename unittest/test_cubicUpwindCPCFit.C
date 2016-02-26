#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "TestPolynomialFit.H"
#include "TestStencils.H"

Approx approx = Approx::custom().epsilon(0.001);

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    const Foam::label faceI = 11;

    Test::PolynomialFit fit(twelvePointStencil(), faceI);
    CHECK(fit.coefficients()[0] == approx(0.875));
    CHECK(fit.coefficients()[1] == approx(0.297));
    CHECK(fit.coefficients()[2] == approx(-0.141));
    CHECK(fit.coefficients()[3] == approx(-0.031));
    CHECK(fit.coefficients()[4] == approx(0.031));
    CHECK(fit.coefficients()[5] == approx(0.008));
    CHECK(fit.coefficients()[6] == approx(-0.086));
    CHECK(fit.coefficients()[7] == approx(0.047));
    CHECK(fit.coefficients()[8] == approx(0.031));
    CHECK(fit.coefficients()[9] == approx(0.008));
    CHECK(fit.coefficients()[10] == approx(-0.086));
    CHECK(fit.coefficients()[11] == approx(0.047));
}

TEST_CASE("fit full-size stencil to uniform set of points in local coords")
{
    Test::PolynomialFit fit(twelvePointStencil());
    CHECK(fit.coefficients()[0] == approx(0.875));
    CHECK(fit.coefficients()[1] == approx(0.297));
    CHECK(fit.coefficients()[2] == approx(-0.141));
    CHECK(fit.coefficients()[3] == approx(-0.031));
    CHECK(fit.coefficients()[4] == approx(0.031));
    CHECK(fit.coefficients()[5] == approx(0.008));
    CHECK(fit.coefficients()[6] == approx(-0.086));
    CHECK(fit.coefficients()[7] == approx(0.047));
    CHECK(fit.coefficients()[8] == approx(0.031));
    CHECK(fit.coefficients()[9] == approx(0.008));
    CHECK(fit.coefficients()[10] == approx(-0.086));
    CHECK(fit.coefficients()[11] == approx(0.047));
}

// while cubicUpwindCPCFit is a correction on upwind, 
// I've added a test for linear correction to help ensure
// that we don't break those particular code paths
TEST_CASE("fit using linear correction")
{

    bool linearCorrection = true;
    Test::PolynomialFit fit(twelvePointStencil(), linearCorrection);
    CHECK(fit.coefficients()[0] == approx(0.875));
    CHECK(fit.coefficients()[1] == approx(0.297));
    CHECK(fit.coefficients()[2] == approx(-0.141));
    CHECK(fit.coefficients()[3] == approx(-0.031));
    CHECK(fit.coefficients()[4] == approx(0.031));
    CHECK(fit.coefficients()[5] == approx(0.008));
    CHECK(fit.coefficients()[6] == approx(-0.086));
    CHECK(fit.coefficients()[7] == approx(0.047));
    CHECK(fit.coefficients()[8] == approx(0.031));
    CHECK(fit.coefficients()[9] == approx(0.008));
    CHECK(fit.coefficients()[10] == approx(-0.086));
    CHECK(fit.coefficients()[11] == approx(0.047));
}

//TEST_CASE("dimensions of FixedPolynomialMatrix")
//{
//}
