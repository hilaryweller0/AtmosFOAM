#include "WeightedMatrix.H"
#include "SVD.H"

Foam::WeightedMatrix::WeightedMatrix(const scalarRectangularMatrix B) : B(B) {};

scalarRectangularMatrix Foam::WeightedMatrix::pseudoInverse() const
{
    SVD svd(B, SMALL);
    return svd.VSinvUt();
}

void Foam::WeightedMatrix::applyStencilPointWeights(const scalarList& weights)
{
    for (label i = 0; i < B.n(); i++)
    {
        for (label j = 0; j < B.m(); j++)
        {
            B[i][j] *= weights[i];
        }
    }
}

void Foam::WeightedMatrix::multiplyConstantAndLinearWeights(const scalar weight)
{
    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= weight;
        B[i][1] *= weight;
    }
}

void Foam::WeightedMatrix::multiplyUpwindWeight(const scalar weight)
{
    multiplyStencilPointWeight(0, weight);
}

void Foam::WeightedMatrix::multiplyDownwindWeight(const scalar weight)
{
    multiplyStencilPointWeight(1, weight);
}

void Foam::WeightedMatrix::multiplyStencilPointWeight(const label index, const scalar weight)
{
    for (label j = 0; j < B.m(); j++)
    {
        B[index][j] *= weight;
    }
}
