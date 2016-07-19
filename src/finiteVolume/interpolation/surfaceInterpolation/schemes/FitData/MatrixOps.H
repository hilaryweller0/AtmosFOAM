#ifndef MatrixOps_H
#define MatrixOps_H

#include <limits>
#include "fvCFD.H"

#include <Eigen/Core>

namespace Foam
{
    template<typename Matrix_Type>
    class MatrixOps
    {
    public:
        static Matrix_Type pseudoInverse(const Matrix_Type &a, double epsilon = std::numeric_limits<double>::epsilon());
    };
}

#ifdef NoRepository
#   include "MatrixOps.C"
#endif

#endif
