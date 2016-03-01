#ifndef PolynomialMatrix_H
#define PolynomialMatrix_H

namespace Foam
{
    template<class Polynomial>
    class PolynomialMatrix
    {
    public:
        virtual ~PolynomialMatrix() {};

        scalarRectangularMatrix pseudoInverse() const;

        virtual scalarRectangularMatrix matrix() const = 0;
    };
}

#ifdef NoRepository
#   include "PolynomialMatrix.C"
#endif

#endif
