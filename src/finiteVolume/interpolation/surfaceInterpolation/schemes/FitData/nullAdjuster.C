#include "nullAdjuster.H"

Foam::nullAdjuster::nullAdjuster
(
    polynomialMatrix& matrix,
    scalarList& coefficients,
    scalarList& wts
)
:
    matrix(matrix),
    coefficients(coefficients),
    wts(wts)
{}

bool Foam::nullAdjuster::adjustWeights()
{
    scalarRectangularMatrix Binv = matrix.pseudoInverse();

    for (label i=0; i<coefficients.size(); i++)
    {
        coefficients[i] = wts[0]*wts[i]*Binv[0][i];
    }

    return true;
}
