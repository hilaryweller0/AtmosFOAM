#include "nullAdjuster.H"

Foam::nullAdjuster::nullAdjuster
(
    weightedMatrix& matrix,
    scalarList& coefficients,
    fitWeights& weights
)
:
    matrix(matrix),
    coefficients(coefficients),
    weights(weights)
{}

bool Foam::nullAdjuster::adjustWeights()
{
    scalarRectangularMatrix Binv = matrix.pseudoInverse();

    for (label i=0; i<coefficients.size(); i++)
    {
        coefficients[i] = weights.constant()*weights[i]*Binv[0][i];
    }

    return true;
}
