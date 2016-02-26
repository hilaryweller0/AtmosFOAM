#include "FixedPolynomialMatrix.H"

template<class Polynomial>
Foam::FixedPolynomialMatrix<Polynomial>::FixedPolynomialMatrix(
        const List<point>& stencil,
        const direction dimensions)
:
    B(stencil.size(), Polynomial::nTerms(dimensions), scalar(0)),
    dimensions(dimensions)
{
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::setStencilPoint(
        const label index, const point& p)
{
        Polynomial::addCoeffs(B[index], p, 1, dimensions);
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::applyStencilPointWeights(const scalarList& weights)
{
    for (label i = 0; i < B.n(); i++)
    {
        for (label j = 0; j < B.m(); j++)
        {
            B[i][j] *= weights[i];
        }
    }
}
