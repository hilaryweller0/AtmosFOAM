#ifndef unweightedMatrix_H
#define unweightedMatrix_H

#include "fvCFD.H"
#include "polynomialMatrix.H"

namespace Foam
{
    class unweightedMatrix : public polynomialMatrix
    {
    public:
        unweightedMatrix(scalarRectangularMatrix& B);

        void applyStencilPointWeights(const scalarList& weights);
        void multiplyConstantAndLinearWeights(const scalar weight);
        void multiplyUpwindWeight(const scalar weight);
        void multiplyDownwindWeight(const scalar weight);
        scalarRectangularMatrix pseudoInverse() const;
    };
}

#endif

