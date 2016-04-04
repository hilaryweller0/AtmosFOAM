#include "stabiliser.H"

// TODO: maybe change signature to accept the fitWeights and scalarRectangularMatrix
// then it will be easier to construct weightedMatrices ourselves
bool Foam::stabiliser::stabilise
(
    const weightedMatrix& matrix,
    fitCoefficients& coefficients
) const
{
    fitCoefficients c(coefficients);

    matrix.populate(c);

//    Info << "*** coeffs (" << matrix.columns() << " terms) " << c << endl;

    label columns = matrix.columns() - 1;

    while (!c.stable() && columns > 0)
    {
        autoPtr<weightedMatrix> m = matrix.truncateTo(columns);
        m->populate(c);
//        Info << "*** coeffs (" << columns << " terms) " << c << endl;
        columns--;
    }

    if (!c.stable())
    {
        // create new weightedMatrix with upwind weight, no downwind weight
        // eventually, we could try removing the downwind weight gradually (from 5 to 1)
        // stopping when stability is achieved.  this should give better accuracy
        // because we've seen in test results that errors are smaller when central weights
        // are used compared to just upwind weighting.
    }

    coefficients.copyFrom(c);

    return c.stable();
}
