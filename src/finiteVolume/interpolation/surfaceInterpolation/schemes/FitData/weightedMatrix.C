#include "weightedMatrix.H"

#include "MatrixOps.H"

#include <Eigen/Core>

using Eigen::MatrixXd;

Foam::weightedMatrix::weightedMatrix
(
    const scalarRectangularMatrix& matrix
)
:
    B(matrix.m(), matrix.n()),
    weights(new fitWeights(matrix.m()))
{
    this->B = matrix;
}

Foam::weightedMatrix::weightedMatrix
(
    const scalarRectangularMatrix& matrix,
    const fitWeights& weights,
    bool applyWeights
)
:
    B(matrix.m(), matrix.n()),
    weights(new fitWeights(weights))
{
    this->B = matrix;

    if (applyWeights)
    {
        applyStencilPointWeights(this->weights);

        for (label i = 0; i < this->B.n(); i++)
        {
            this->B[i][0] *= this->weights->constant();
            this->B[i][1] *= this->weights->xLinear();
        }
    }
}

Foam::weightedMatrix::weightedMatrix
(
    const weightedMatrix& matrix,
    const fitWeights& weights
)
:
    B(matrix.B.m(), matrix.B.n()),
    weights(new fitWeights(weights))
{
    this->B = matrix.B;

    applyStencilPointWeights(this->weights);

    for (label i = 0; i < this->B.n(); i++)
    {
        this->B[i][0] *= this->weights->constant();
        this->B[i][1] *= this->weights->xLinear();
    }
}

void Foam::weightedMatrix::populate(fitCoefficients& coefficients) const
{
	MatrixXd m(B.n(), B.m());
    for (label i = 0; i < B.n(); i++)
    {
        for (label j = 0; j < B.m(); j++)
        {
            m(i,j) = B[i][j];
        }
    }

    MatrixXd Binv = MatrixOps<MatrixXd>::pseudoInverse(m);

    for (label i=0; i<coefficients.size(); i++)
    {
        coefficients[i] = weights->constant()*weights()[i]*Binv(0,i);
    }
}

autoPtr<weightedMatrix> Foam::weightedMatrix::subset(const labelList columns) const
{
    scalarRectangularMatrix subset(B.n(), columns.size());

    for (label i = 0; i < subset.n(); i++)
    {
        for (label j = 0; j < subset.m(); j++)
        {
            subset[i][j] = B[i][columns[j]];
        }
    }

    return autoPtr<weightedMatrix>
    (
        new weightedMatrix(subset, weights, false)
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

