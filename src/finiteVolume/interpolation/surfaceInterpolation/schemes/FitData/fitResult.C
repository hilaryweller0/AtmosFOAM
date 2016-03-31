#include "fitResult.H"

Foam::fitResult::fitResult(
        const fitCoefficients& coefficients,
        const bool good,
        const label polynomialTerms
)
:
    coefficients(coefficients),
    good(good),
    polynomialTerms(polynomialTerms)
{};
