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
    expectedCoefficients[ 0] =  0.5804;
    expectedCoefficients[ 1] =  0.4198;
    expectedCoefficients[ 2] = -0.0513;
    expectedCoefficients[ 3] = -0.0531;
    expectedCoefficients[ 4] =  0.0486;
    expectedCoefficients[ 5] =  0.0556;

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
    expectedCoefficients[ 0] =  0.9965459942;
    expectedCoefficients[ 1] =  0.1344156036;
    expectedCoefficients[ 2] = -0.1449563476;
    expectedCoefficients[ 3] = -0.0004502506;
    expectedCoefficients[ 4] =  0.0144450004;

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

TEST_CASE("eightPointsWithDiagonal")
{
    List<point> stencilPoints(8, point(0, 0, 0));
    stencilPoints[0] = point(-1, 0.0425922049484, 0);
    stencilPoints[1] = point(1.05781693914, -0.0369040044104, 0);
    stencilPoints[2] = point(-5.10831948925, 0.22706220493, 0);
    stencilPoints[3] = point(-3.05606764231, 0.133444521133, 0);
    stencilPoints[4] = point(-5.17064352149, 1.10939830377, 0);
    stencilPoints[5] = point(-3.1023861129, 1.10939830377, 0);
    stencilPoints[6] = point(-1.0341287043, 1.10939830377, 0);
    stencilPoints[7] = point(1.0341287043, 1.10939830377, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.8507;
    expectedCoefficients[ 1] =  0.3320;
    expectedCoefficients[ 2] =  0.0227;
    expectedCoefficients[ 3] = -0.2036;
    expectedCoefficients[ 4] =  0.0385;
    expectedCoefficients[ 5] = -0.1051;
    expectedCoefficients[ 6] =  0.0971;
    expectedCoefficients[ 7] = -0.0323;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("eightPointsWithDiagonal2")
{
    List<point> stencilPoints(8, point(0, 0, 0));
    stencilPoints[0] = point(-1, -0.0669388588762, 0);
    stencilPoints[1] = point(0.833644437615, 0.0794066407836, 0);
    stencilPoints[2] = point(-4.75714056459, -0.217550549489, 0);
    stencilPoints[3] = point(-2.8731750957, -0.165967824541, 0);
    stencilPoints[4] = point(-4.74883733854, 0.906972333924, 0);
    stencilPoints[5] = point(-2.84930240312, 0.906972333924, 0);
    stencilPoints[6] = point(-0.949767467708, 0.906972333924, 0);
    stencilPoints[7] = point(0.949767467708, 0.906972333924, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.7690;
    expectedCoefficients[ 1] =  0.4155;
    expectedCoefficients[ 2] =  0.0338;
    expectedCoefficients[ 3] = -0.2085;
    expectedCoefficients[ 4] =  0.0250;
    expectedCoefficients[ 5] = -0.0878;
    expectedCoefficients[ 6] =  0.1040;
    expectedCoefficients[ 7] = -0.0510;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("twoByThreeVerySmallDiagonal")
{
    List<point> stencilPoints(6, point(0, 0, 0));
    stencilPoints[0] = point(-1, -0.000810356386742, 0);
    stencilPoints[1] = point(1.00145581289, -8.74386005651e-14, 0);
    stencilPoints[2] = point(-0.993588838779, 3.33540741754, 0);
    stencilPoints[3] = point(-1.00145581289, -1.66909302149, 0);
    stencilPoints[4] = point(1.00145581289, -1.66909302149, 0);
    stencilPoints[5] = point(1.00145581289, 3.33818604298, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.5005;
    expectedCoefficients[ 1] =  0.4996;
    expectedCoefficients[ 2] =  0.0000;
    expectedCoefficients[ 3] = -0.0002;
    expectedCoefficients[ 4] =  0.0000;
    expectedCoefficients[ 5] =  0.0000;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}

TEST_CASE("twoByThreeDiagonal")
{
    List<point> stencilPoints(6, point(0, 0, 0));
    stencilPoints[0] = point(-1, 0.0149698735917, 0);
    stencilPoints[1] = point(1.02740059135, -2.24259694664e-14, 0);
    stencilPoints[2] = point(-0.907001583868, -3.3816285968, 0);
    stencilPoints[3] = point(-1.0162149057, 3.41915929803, 0);
    stencilPoints[4] = point(1.02740059135, -3.42466863784, 0);
    stencilPoints[5] = point(1.02740059135, 3.42466863784, 0);

    const localStencil stencil(stencilPoints);
    fitCoefficients actualCoefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    fitCoefficients expectedCoefficients(stencil.size(), false, 0);
    expectedCoefficients[ 0] =  0.5067;
    expectedCoefficients[ 1] =  0.4932;
    expectedCoefficients[ 2] =  0.0012;
    expectedCoefficients[ 3] = -0.0011;
    expectedCoefficients[ 4] = -0.0001;
    expectedCoefficients[ 5] =  0.0000;

    autoPtr<fitResult> actual = fitPolynomial(actualCoefficients, weights, stencil);

    check(actualCoefficients, expectedCoefficients);
}
