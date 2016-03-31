#include "stabiliser.H"

bool Foam::stabiliser::stabilise
(
    const weightedMatrix& matrix,
    fitCoefficients& coefficients
) const
{
    fitCoefficients c(coefficients);

    matrix.populate(c);

    if (!c.stable())
    {
        // remove degree-3 columns from matrix
        autoPtr<weightedMatrix> m = matrix.truncateToAtMost(6); // columns
        m->populate(c);
    }

    coefficients.copyFrom(c);

    return true;
}
