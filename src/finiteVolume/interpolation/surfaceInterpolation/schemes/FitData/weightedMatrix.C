#include "weightedMatrix.H"
#include "SVD.H"

Foam::weightedMatrix::weightedMatrix(scalarRectangularMatrix& B) : B(B) {};

void Foam::weightedMatrix::apply(const fitWeights& weights)
{
    applyStencilPointWeights(weights);

    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= weights.constant();
        B[i][1] *= weights.xLinear();
    }
}

scalarRectangularMatrix Foam::weightedMatrix::pseudoInverse() const
{
    SVD svd(B, SMALL);
    return svd.VSinvUt();
}

void Foam::weightedMatrix::applyStencilPointWeights(const fitWeights& weights)
{
    for (label i = 0; i < B.n(); i++)
    {
        for (label j = 0; j < B.m(); j++)
        {
            B[i][j] *= weights[i];
        }
    }
}

