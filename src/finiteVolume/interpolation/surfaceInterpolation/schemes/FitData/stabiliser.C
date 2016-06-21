#include <vector>

#include "stabiliser.H"
#include "weightedMatrix.H"

label Foam::stabiliser::stabilise
(
    const scalarRectangularMatrix& B,
    fitWeights& weights,
    fitCoefficients& coefficients,
    label faceI,
    bool owner
) const
{
    autoPtr<weightedMatrix> unweightedMatrix = findStabilisableMatrix(B, weights, coefficients, faceI, owner);
    fitCoefficients c(coefficients);

    label columns;

    do 
    {
        weightedMatrix matrix(unweightedMatrix, weights);
        columns = matrix.columns();
        matrix.populate(c);

        weights.downwind() -= 1.0;
    } while (!c.stable() && weights.downwind() >= 1.0);

    weights.downwind() += 1.0;

    // revert to something akin to linearUpwind
    if (!c.stable())
    {
        weights.downwind() = 1.0;
        weightedMatrix matrix(B, weights);
        // include constant and linear terms
        labelList columnIndices(0, 0);
        columnIndices.append(0);
        columnIndices.append(1);
        columnIndices.append(2); // FIXME: hardwired for 2D stencils
        columns = 3;
        autoPtr<weightedMatrix> m = matrix.subset(columnIndices);
        m->populate(c);
    }

    coefficients.copyFrom(c);

    return c.stable() ? columns : 0;
}

autoPtr<weightedMatrix> Foam::stabiliser::findStabilisableMatrix
(
    const scalarRectangularMatrix& B,
    const fitWeights& weights,
    const fitCoefficients& coefficients,
    label faceI,
    bool owner
) const
{
    fitCoefficients c(coefficients);
    weightedMatrix matrix(B, weights);

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
            if (c.central_are_largest())
            {
                Info << "*** polyTerms " << columnIndices << endl;
                weightedMatrix unweighted(B);
                return unweighted.subset(columnIndices);
            }
        } while (std::prev_permutation(v.begin(), v.end()));
    }

    // FIXME: oh bugger! raise exception?
    Warning << "Well crap!" << endl;
    return autoPtr<weightedMatrix>(new weightedMatrix(B));
}
