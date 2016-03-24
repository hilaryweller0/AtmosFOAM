#ifndef polynomialMatrix_H
#define polynomialMatrix_H

#include "fvCFD.H"

namespace Foam
{
    class polynomialMatrix
    {
    public:
        polynomialMatrix(scalarRectangularMatrix& B);
        virtual ~polynomialMatrix() {};

        virtual void applyStencilPointWeights(const scalarList& weights) = 0;
        virtual void multiplyConstantAndLinearWeights(const scalar weight) = 0;
        virtual void multiplyUpwindWeight(const scalar weight) = 0;
        virtual void multiplyDownwindWeight(const scalar weight) = 0;
        virtual scalarRectangularMatrix pseudoInverse() const = 0;

    protected:
        scalarRectangularMatrix& B;
    };
}

#endif
