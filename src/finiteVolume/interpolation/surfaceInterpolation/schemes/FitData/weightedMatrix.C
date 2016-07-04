#include "weightedMatrix.H"

#include <Eigen/SVD>

using Eigen::MatrixXd;
using Eigen::JacobiSVD;

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
//
// from http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257#c14
template<typename _Matrix_Type_>
_Matrix_Type_ Foam::weightedMatrix::pseudoInverse(const _Matrix_Type_ &a, double epsilon) const
{
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
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

    MatrixXd Binv = pseudoInverse(m);

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

