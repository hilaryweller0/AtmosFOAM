#include "testablePolynomialFit.H"
#include "PolynomialFit2.H"
#include "cubicUpwindCPCFitPolynomial.H"

autoPtr<fitResult> fitPolynomial
(
        fitCoefficients& coefficients,
        fitWeights& weights,
        const localStencil& stencil
)
{
    const direction dimensions = 2;
    PolynomialFit2<cubicUpwindCPCFitPolynomial> polynomialFit(dimensions);
    autoPtr<fitResult> result = polynomialFit.fit(coefficients, weights, stencil);
    coefficients[0] += 1;
    return result;
}
