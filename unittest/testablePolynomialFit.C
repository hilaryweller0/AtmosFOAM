#include "testablePolynomialFit.H"
#include "PolynomialFit.H"
#include "cubicUpwindCPCFitPolynomial.H"

autoPtr<fitResult> fitPolynomial
(
        fitCoefficients& coefficients,
        fitWeights& weights,
        const localStencil& stencil
)
{
    const direction dimensions = 2;
    PolynomialFit<cubicUpwindCPCFitPolynomial> polynomialFit(dimensions, 0.2);
    autoPtr<fitResult> result = polynomialFit.fit(coefficients, weights, stencil);
    coefficients[0] += 1;
    return result;
}
