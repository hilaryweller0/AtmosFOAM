#include "MatrixOps.H"

#include <Eigen/SVD>

using Eigen::MatrixXd;
using Eigen::JacobiSVD;

// from http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257#c14
template<typename Matrix_Type>
Matrix_Type Foam::MatrixOps<Matrix_Type>::pseudoInverse(const Matrix_Type &a, double epsilon)
{
/*    Info << a.rows() << "x" << a.cols() << endl;
    for (label row=0; row<a.rows(); row++)
    {
        for (label col=0; col<a.cols(); col++)
        {
            Info << a(row, col) << " ";
        }
        Info << endl;
    }
    */
    JacobiSVD<Matrix_Type> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    Matrix_Type inv = svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
/*
    Info << inv.rows() << "x" << inv.cols() << endl;
    for (label row=0; row<inv.rows(); row++)
    {
        for (label col=0; col<inv.cols(); col++)
        {
            Info << inv(row, col) << " ";
        }
        Info << endl;
    }
*/
    return inv;
}
