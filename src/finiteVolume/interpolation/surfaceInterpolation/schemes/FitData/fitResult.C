#include "fitResult.H"

Foam::fitResult::fitResult(
        const List<point>& stencilPoints,
        const fitCoefficients& coefficients,
        const bool good,
        const label polynomialTerms
)
:
    stencilPoints(stencilPoints),
    coefficients(coefficients),
    good(good),
    polynomialTerms(polynomialTerms)
{};
