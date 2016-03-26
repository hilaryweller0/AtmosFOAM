#include "weightedMatrix.H"
#include "SVD.H"

Foam::weightedMatrix::weightedMatrix(scalarRectangularMatrix& B) : B(B) {};

scalarRectangularMatrix Foam::weightedMatrix::pseudoInverse() const
{
    SVD svd(B, SMALL);
    return svd.VSinvUt();
}

void Foam::weightedMatrix::applyStencilPointWeights(const scalarList& weights)
{
    for (label i = 0; i < B.n(); i++)
    {
        for (label j = 0; j < B.m(); j++)
        {
            B[i][j] *= weights[i];
        }
    }
}

void Foam::weightedMatrix::multiplyConstantAndLinearWeights(const scalar weight)
{
    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= weight;
        B[i][1] *= weight;
    }
}

void Foam::weightedMatrix::multiplyUpwindWeight(const scalar weight)
{
    multiplyStencilPointWeight(0, weight);
}

void Foam::weightedMatrix::multiplyDownwindWeight(const scalar weight)
{
    multiplyStencilPointWeight(1, weight);
}

void Foam::weightedMatrix::multiplyStencilPointWeight(const label index, const scalar weight)
{
    for (label j = 0; j < B.m(); j++)
    {
        B[index][j] *= weight;
    }
}
