#include "nullAdjuster.H"

Foam::nullAdjuster::nullAdjuster
(
    weightedMatrix& matrix,
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
        // wts[0] is for the constant term weighting
        // we should delegate to obtain this value

        // wts[i] is for the stencil cell weightings
        // and we should similarly delegate to obtain these values
        coefficients[i] = wts[0]*wts[i]*Binv[0][i];
    }

    return true;
}
