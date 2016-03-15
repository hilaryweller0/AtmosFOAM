#include "Fit.H"

Foam::Fit::Fit(
        const List<point>& stencilPoints,
        const scalarList& coeffs,
        const bool good
)
:
    stencilPoints(stencilPoints),
    coeffs(coeffs),
    good(good)
{};
