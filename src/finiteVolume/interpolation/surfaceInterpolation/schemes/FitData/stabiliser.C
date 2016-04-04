#include "stabiliser.H"
#include "weightedMatrix.H"

// TODO: maybe change signature to accept the fitWeights and scalarRectangularMatrix
// then it will be easier to construct weightedMatrices ourselves
bool Foam::stabiliser::stabilise
(
    const scalarRectangularMatrix& B,
    fitWeights& weights,
    fitCoefficients& coefficients
) const
{
    fitCoefficients c(coefficients);

    weightedMatrix matrix(B, weights);
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

    while (!c.stable() && weights.downwind() >= 2.0)
    {
        weights.downwind() -= 1.0;
        weightedMatrix matrix(B, weights);
        matrix.populate(c);
        Info << "*** coeffs (" << weights.downwind() << " downwind weight) " << c << endl;
    }

    coefficients.copyFrom(c);

    return c.stable();
}
