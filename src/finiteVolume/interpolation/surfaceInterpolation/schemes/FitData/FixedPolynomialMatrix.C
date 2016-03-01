#include "FixedPolynomialMatrix.H"

template<class Polynomial>
Foam::FixedPolynomialMatrix<Polynomial>::FixedPolynomialMatrix(
        const List<point>& stencil,
        const direction dimensions)
:
    B(stencil.size(), Polynomial::nTerms(dimensions), scalar(0))
{
    forAll(stencil, i)
    {
        Polynomial::addCoeffs(B[i], stencil[i], 1, dimensions);
    }
}

template<class Polynomial>
scalarRectangularMatrix Foam::FixedPolynomialMatrix<Polynomial>::matrix() const
{
    return B;
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
    multiplyStencilPointWeight(0, weight);
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::multiplyDownwindWeight(const scalar weight)
{
    multiplyStencilPointWeight(1, weight);
}

template<class Polynomial>
void Foam::FixedPolynomialMatrix<Polynomial>::multiplyStencilPointWeight(const label index, const scalar weight)
{
    for (label j = 0; j < B.m(); j++)
    {
        B[index][j] *= weight;
    }
}
