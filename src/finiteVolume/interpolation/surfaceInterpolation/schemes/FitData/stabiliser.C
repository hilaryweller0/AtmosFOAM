#include "stabiliser.H"
#include "weightedMatrix.H"

label Foam::stabiliser::stabilise
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

    label columns = matrix.columns();

    while (!c.stable() && columns > 1)
    {
        columns--;
        autoPtr<weightedMatrix> m = matrix.truncateTo(columns);
        m->populate(c);
//        Info << "*** coeffs (" << columns << " terms) " << c << endl;
    }

    while (!c.stable() && weights.downwind() >= 2.0)
    {
        columns = matrix.columns();
        weights.downwind() -= 1.0;
        weightedMatrix matrix(B, weights);
        matrix.populate(c);
//        Info << "*** coeffs (" << weights.downwind() << " downwind weight) " << c << endl;
    }

    coefficients.copyFrom(c);

    return c.stable() ? columns : 0;
}
