#ifndef PolynomialMatrix_H
#define PolynomialMatrix_H

namespace Foam
{
    template<class Polynomial>
    class AbstractPolynomial
    {
    public:
        virtual ~AbstractPolynomial() {};

        virtual scalarRectangularMatrix matrix() const = 0;
    };
}

#endif
