#include "weightedMatrix.H"
#include "SVD.H"

Foam::weightedMatrix::weightedMatrix
(
    scalarRectangularMatrix& B,
    const fitWeights& weights
)
:
    B(B),
    weights(weights)
{
    applyStencilPointWeights(weights);

    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= weights.constant();
        B[i][1] *= weights.xLinear();
    }
}

void Foam::weightedMatrix::populate(fitCoefficients& coefficients) const
{
    SVD svd(B, SMALL);
    const scalarRectangularMatrix& Binv = svd.VSinvUt();

    for (label i=0; i<coefficients.size(); i++)
    {
        coefficients[i] = weights.constant()*weights[i]*Binv[0][i];
    }
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

