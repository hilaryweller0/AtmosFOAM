#ifndef AbstractPolynomial_H
#define AbstractPolynomial_H

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
