#include "fitResult.H"

Foam::fitResult::fitResult
(
        const localStencil stencil,
        const fitCoefficients& coefficients,
        const fitWeights& weights,
        const bool good,
        const uint32_t polynomial,
        const label polynomialTerms
)
:
    stencil(stencil),
    coefficients(coefficients),
    weights(weights),
    good(good),
    polynomial(polynomial),
    polynomialTerms(polynomialTerms)
{};
