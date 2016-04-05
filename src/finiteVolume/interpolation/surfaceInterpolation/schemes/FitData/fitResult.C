#include "fitResult.H"

Foam::fitResult::fitResult
(
        const localStencil stencil,
        const fitCoefficients& coefficients,
        const fitWeights& weights,
        const bool good,
        const label polynomialTerms
)
:
    stencil(stencil),
    coefficients(coefficients),
    weights(weights),
    good(good),
    polynomialTerms(polynomialTerms)
{};
