#include "PolynomialMatrix.H"
#include "SVD.H"

template<class Polynomial>
scalarRectangularMatrix Foam::PolynomialMatrix<Polynomial>::pseudoInverse() const
{
    SVD svd(matrix(), SMALL);
    return svd.VSinvUt();
}
