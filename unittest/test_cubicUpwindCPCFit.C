#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include <assert.h>
#include "PolynomialFit2.H"
#include "cubicUpwindCPCFitPolynomial.H"
#include "checks.H"

TEST_CASE("uniform2DQuadInterior")
{
    List<point> stencilPoints(12, point(0, 0, 0));
    stencilPoints[0] = point(-1, 0, 0);
    stencilPoints[1] = point(1, 0, 0);
    stencilPoints[2] = point(-3, 0, 0);
    stencilPoints[3] = point(-5, 0, 0);
    stencilPoints[4] = point(-1, -2, 0);
    stencilPoints[5] = point(1, -2, 0);
    stencilPoints[6] = point(-3, -2, 0);
    stencilPoints[7] = point(-5, -2, 0);
    stencilPoints[8] = point(-1, 2, 0);
    stencilPoints[9] = point(1, 2, 0);
    stencilPoints[10] = point(-3, 2, 0);
    stencilPoints[11] = point(-5, 2, 0); 

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    const direction dimensions = 2;
    PolynomialFit2<cubicUpwindCPCFitPolynomial> polynomialFit(dimensions);
    polynomialFit.fit(actualCoefficients, weights, stencil);
    actualCoefficients[0] += 1;

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.8745;
    expectedCoefficients[ 1] =  0.2969;
    expectedCoefficients[ 2] = -0.1406;
    expectedCoefficients[ 3] = -0.0312;
    expectedCoefficients[ 4] =  0.0313;
    expectedCoefficients[ 5] =  0.0078;
    expectedCoefficients[ 6] = -0.0859;
    expectedCoefficients[ 7] =  0.0469;
    expectedCoefficients[ 8] =  0.0313;
    expectedCoefficients[ 9] =  0.0078;
    expectedCoefficients[10] = -0.0859;
    expectedCoefficients[11] =  0.0469;

    check(actualCoefficients, expectedCoefficients);
}
