#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "TestPolynomialFit.H"

Approx approx = Approx::custom().epsilon(0.001);

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    Foam::List<point> twelvePointStencil(12, point(0, 0, 0));
    twelvePointStencil[0] = point(-0.5, 0, 0);
    twelvePointStencil[1] = point(0.5, 0, 0);
    twelvePointStencil[2] = point(-1.5, 0, 0);
    twelvePointStencil[3] = point(-2.5, 0, 0);
    twelvePointStencil[4] = point(-0.5, 0, 1);
    twelvePointStencil[5] = point(0.5, 0, 1);
    twelvePointStencil[6] = point(-1.5, 0, 1);
    twelvePointStencil[7] = point(-2.5, 0, 1);
    twelvePointStencil[8] = point(-0.5, 0, -1);
    twelvePointStencil[9] = point(0.5, 0, -1);
    twelvePointStencil[10] = point(-1.5, 0, -1);
    twelvePointStencil[11] = point(-2.5, 0, -1);

    const Foam::label faceI = 11;

    Test::PolynomialFit fit(twelvePointStencil, faceI);
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
    Foam::List<point> twelvePointStencil(12, point(0, 0, 0));
    twelvePointStencil[0] = point(-0.5, 0, 0);
    twelvePointStencil[1] = point(0.5, 0, 0);
    twelvePointStencil[2] = point(-1.5, 0, 0);
    twelvePointStencil[3] = point(-2.5, 0, 0);
    twelvePointStencil[4] = point(-0.5, 0, 1);
    twelvePointStencil[5] = point(0.5, 0, 1);
    twelvePointStencil[6] = point(-1.5, 0, 1);
    twelvePointStencil[7] = point(-2.5, 0, 1);
    twelvePointStencil[8] = point(-0.5, 0, -1);
    twelvePointStencil[9] = point(0.5, 0, -1);
    twelvePointStencil[10] = point(-1.5, 0, -1);
    twelvePointStencil[11] = point(-2.5, 0, -1);

    Test::PolynomialFit fit(twelvePointStencil);
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
