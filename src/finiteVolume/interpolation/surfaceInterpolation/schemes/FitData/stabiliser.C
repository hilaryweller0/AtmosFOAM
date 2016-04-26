#include <vector>

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

//    Info << "*** starting with " << B.m() << " candidate terms" << endl;

    weightedMatrix matrix(B, weights);

    // try combinations of columns to find the best stabilisable matrix
    for (int termCount = B.m(); termCount > 0; termCount--)
    {
        std::vector<bool> v(B.m());
        std::fill(v.begin(), v.end() - B.m() + termCount, true);

        do
        {
            labelList columnIndices(0, 0);
            for (int i = 0; i < B.m(); ++i)
            {
                if (v[i])
                {
                    columnIndices.append(i);
                }
            }
            autoPtr<weightedMatrix> m = matrix.subset(columnIndices);
            m->populate(c);
            //Info << "columns " << columnIndices << " is " << (c.central_are_largest() ? "stabilisable" : "not stabilisable") << endl;
        } while (std::prev_permutation(v.begin(), v.end()));
        Info << endl;
    }

    matrix.populate(c);
    label columns = matrix.columns();

    while (!c.stable() && columns > 1)
    {
        columns--;
        autoPtr<weightedMatrix> m = matrix.truncateTo(columns);
        m->populate(c);
    }

    while (!c.stable() && weights.downwind() >= 2.0)
    {
        columns = matrix.columns();
        weights.downwind() -= 1.0;
        weightedMatrix matrix(B, weights);
        matrix.populate(c);
    }

    coefficients.copyFrom(c);

    return c.stable() ? columns : 0;
}
