#include "unweightedMatrix.H"
#include "SVD.H"

Foam::unweightedMatrix::unweightedMatrix(const scalarRectangularMatrix B) : polynomialMatrix(B) {};

scalarRectangularMatrix Foam::unweightedMatrix::pseudoInverse() const
{
    SVD svd(B, SMALL);
    return svd.VSinvUt();
}

void Foam::unweightedMatrix::applyStencilPointWeights(const scalarList& weights)
{
}

void Foam::unweightedMatrix::multiplyConstantAndLinearWeights(const scalar weight)
{
}

void Foam::unweightedMatrix::multiplyUpwindWeight(const scalar weight)
{
}

void Foam::unweightedMatrix::multiplyDownwindWeight(const scalar weight)
{
}
