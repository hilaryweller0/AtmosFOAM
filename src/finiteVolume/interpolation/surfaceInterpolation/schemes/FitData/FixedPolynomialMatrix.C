#include "FixedPolynomialMatrix.H"
#include "SVD.H"

template<class Polynomial>
Foam::FixedPolynomialMatrix<Polynomial>::FixedPolynomialMatrix(
        const List<point>& stencil,
        const direction dimensions)
:
    B(stencil.size(), Polynomial::nTerms(dimensions), scalar(0)),
    dimensions(dimensions)
{
    forAll(stencil, i)
    {
        Polynomial::addCoeffs(B[i], stencil[i], 1, dimensions);
    }
}

template<class Polynomial>
scalarRectangularMatrix Foam::FixedPolynomialMatrix<Polynomial>::pseudoInverse()
{
    SVD svd(B, SMALL);
    return svd.VSinvUt();
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

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::multiplyConstantAndLinearWeights(const scalar weight)
{
    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= weight;
        B[i][1] *= weight;
    }
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::multiplyUpwindWeight(const scalar weight)
{
    for (label j = 0; j < B.m(); j++)
    {
        B[0][j] *= weight;
    }
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::multiplyDownwindWeight(const scalar weight)
{
    for (label j = 0; j < B.m(); j++)
    {
        B[1][j] *= weight;
    }
}
