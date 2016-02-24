#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "PolynomialFit.H"

Approx approx = Approx::custom().epsilon(0.001);

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    Foam::List<point> stencilPoints(12, point(0, 0, 0));
    stencilPoints[0] = point(-0.5, 0, 0);
    stencilPoints[1] = point(0.5, 0, 0);
    stencilPoints[2] = point(-1.5, 0, 0);
    stencilPoints[3] = point(-2.5, 0, 0);
    stencilPoints[4] = point(-0.5, 0, 1);
    stencilPoints[5] = point(0.5, 0, 1);
    stencilPoints[6] = point(-1.5, 0, 1);
    stencilPoints[7] = point(-2.5, 0, 1);
    stencilPoints[8] = point(-0.5, 0, -1);
    stencilPoints[9] = point(0.5, 0, -1);
    stencilPoints[10] = point(-1.5, 0, -1);
    stencilPoints[11] = point(-2.5, 0, -1);

    const Foam::label faceI = 11;

    Test::PolynomialFit fit(stencilPoints, faceI);
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
