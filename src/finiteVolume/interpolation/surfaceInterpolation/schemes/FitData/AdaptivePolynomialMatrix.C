#include "AdaptivePolynomialMatrix.H"

template<class Polynomial>
Foam::AdaptivePolynomialMatrix<Polynomial>::AdaptivePolynomialMatrix(
        const List<point>& stencil,
        const direction dimensions)
:
    stencil(stencil),
    dimensions(dimensions),
    maxTerms(Polynomial::nTerms(dimensions))
{}

template<class Polynomial>
scalarRectangularMatrix Foam::AdaptivePolynomialMatrix<Polynomial>::matrix()
{
    Foam::scalarRectangularMatrix B(stencil.size(), 2, scalar(0));

    for (int i=0; i<B.n(); i++)
    {
        scalar coefficients[maxTerms]; 
        Polynomial::addCoeffs(coefficients, stencil[i], 1, dimensions);
        for (int j=0; j<B.m(); j++)
        {
            B[i][j] = coefficients[j];
        }
    }

    return B;
}
