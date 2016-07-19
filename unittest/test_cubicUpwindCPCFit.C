#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "PolynomialFit2.H"
#include "cubicUpwindCPCFitPolynomial.H"
#include <assert.h>

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
    fitCoefficients coefficients(stencil.size(), false, 0);
    fitWeights weights(stencil.size());

    const direction dimensions = 2;
    PolynomialFit2<cubicUpwindCPCFitPolynomial> polynomialFit(dimensions);
    polynomialFit.fit(coefficients, weights, stencil);
    coefficients[0] += 1;

    Info << coefficients << endl;
    Info << weights << endl;
    // TODO: assert something
}
