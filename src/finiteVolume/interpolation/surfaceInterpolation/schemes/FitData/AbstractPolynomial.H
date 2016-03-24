#ifndef AbstractPolynomial_H
#define AbstractPolynomial_H

#include "autoPtr.H"

namespace Foam
{
    template<class Polynomial>
    class AbstractPolynomial
    {
    public:
        virtual ~AbstractPolynomial() {};

        virtual autoPtr<scalarRectangularMatrix> matrix() const = 0;
    };
}

#endif
