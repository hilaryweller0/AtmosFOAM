#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include <assert.h>
#include "testablePolynomialFit.H"
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

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("sixPointsWithDiagonal")
{
    List<point> stencilPoints(6, point(0, 0, 0));
    stencilPoints[0] = point(-1, 0.0330314, 0);
    stencilPoints[1] = point(0.907926, -2.64e-14, 0);
    stencilPoints[2] = point(-1.13441, 3.05716, 0);
    stencilPoints[3] = point(-0.87524, -2.99221, 0);
    stencilPoints[4] = point(0.907926, 3.02642, 0);
    stencilPoints[5] = point(0.907926, -3.02642, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.5069;
    expectedCoefficients[ 1] =  0.4931;
    expectedCoefficients[ 2] = -0.0169;
    expectedCoefficients[ 3] = -0.0138;
    expectedCoefficients[ 4] =  0.0143;
    expectedCoefficients[ 5] =  0.0164;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("ninePointsWithDiagonal")
{
    List<point> stencilPoints(9, point(0, 0, 0));
    stencilPoints[0] = point(-1, 1.36424e-16, 0);
    stencilPoints[1] = point(1, 5.80712e-14, 0);
    stencilPoints[2] = point(1, 3.33333, 0);
    stencilPoints[3] = point(1, -3.33333, 0);
    stencilPoints[4] = point(-1, 3.33333, 0);
    stencilPoints[5] = point(-1, -3.33333, 0);
    stencilPoints[6] = point(-3.12554, 3.35946, 0);
    stencilPoints[7] = point(-3.02091, 0.0282014, 0);
    stencilPoints[8] = point(-2.91995, -3.30361, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.6864;
    expectedCoefficients[ 1] =  0.4063;
    expectedCoefficients[ 2] = -0.0150;
    expectedCoefficients[ 3] = -0.0157;
    expectedCoefficients[ 4] =  0.0291;
    expectedCoefficients[ 5] =  0.0322;
    expectedCoefficients[ 6] = -0.0136;
    expectedCoefficients[ 7] = -0.0927;
    expectedCoefficients[ 8] = -0.0170;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("threeDownwindTwoUpwind")
{
    List<point> stencilPoints(5, point(0, 0, 0));
    stencilPoints[0] = point(-0.318132398309, 1, 0);
    stencilPoints[1] = point(1.62494112442, -0.0619155571362, 0);
    stencilPoints[2] = point(-0.519892365614, 6.24311965797, 0);
    stencilPoints[3] = point(1.97959459194, 5.68059436497, 0);
    stencilPoints[4] = point(1.6715638405, -5.58574321789, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.8697;
    expectedCoefficients[ 1] =  0.1475;
    expectedCoefficients[ 2] = -0.0366;
    expectedCoefficients[ 3] = -0.0465;
    expectedCoefficients[ 4] =  0.0659;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("threeDownwindTwoUpwind2")
{
    List<point> stencilPoints(5, point(0, 0, 0));
    stencilPoints[0] = point(-0.375692991674, 1, 0);
    stencilPoints[1] = point(2.04670349731, 0.151415600652, 0);
    stencilPoints[2] = point(2.68167823645, -5.90562193016, 0);
    stencilPoints[3] = point(1.48171410264, 6.45067028153, 0);
    stencilPoints[4] = point(-1.05780329301, 6.2993205546, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.9205;
    expectedCoefficients[ 1] =  0.0921;
    expectedCoefficients[ 2] =  0.0695;
    expectedCoefficients[ 3] = -0.0457;
    expectedCoefficients[ 4] = -0.0364;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}
