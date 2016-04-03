#include "fitResult.H"

Foam::fitResult::fitResult(
        const localStencil stencil,
        const fitCoefficients& coefficients,
        const bool good,
        const label polynomialTerms
)
:
    stencil(stencil),
    coefficients(coefficients),
    good(good),
    polynomialTerms(polynomialTerms)
{};
