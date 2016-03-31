#include "weightedMatrix.H"
#include "SVD.H"

Foam::weightedMatrix::weightedMatrix
(
    const scalarRectangularMatrix& matrix,
    const fitWeights& weights,
    bool applyWeights
)
:
    B(matrix.m(), matrix.n()),
    weights(weights)
{
    this->B = matrix;

    if (applyWeights)
    {
        applyStencilPointWeights(weights);

        for (label i = 0; i < this->B.n(); i++)
        {
            this->B[i][0] *= weights.constant();
            this->B[i][1] *= weights.xLinear();
        }
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

autoPtr<weightedMatrix> Foam::weightedMatrix::truncateToAtMost(const label columns) const
{
    scalarRectangularMatrix truncated(B.n(), min(columns, B.m()));

    for (label i = 0; i < truncated.n(); i++)
    {
        for (label j = 0; j < truncated.m(); j++)
        {
            truncated[i][j] = B[i][j];
        }
    }

    return autoPtr<weightedMatrix>
    (
        new weightedMatrix(truncated, weights, false)
    );
}

label Foam::weightedMatrix::columns() const
{
    return B.m();
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

