#include "stabiliser.H"

bool Foam::stabiliser::stabilise
(
    const weightedMatrix& matrix,
    fitCoefficients& coefficients
) const
{
    matrix.populate(coefficients);

    return true;
}
