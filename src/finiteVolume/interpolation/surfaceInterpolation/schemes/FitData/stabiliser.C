#include "stabiliser.H"

bool Foam::stabiliser::stabilise
(
    const weightedMatrix& matrix,
    fitCoefficients& coefficients
) const
{
    fitCoefficients c(coefficients);

    matrix.populate(c);

    label columns = matrix.columns() - 1;

    while (!c.stable() && columns > 0)
    {
        autoPtr<weightedMatrix> m = matrix.truncateToAtMost(columns);
        m->populate(c);
        columns--;
    }

    if (!c.stable()) Info << c << endl;

    coefficients.copyFrom(c);

    return c.stable();
}
