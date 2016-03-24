#ifndef weightedMatrix_H
#define weightedMatrix_H

#include "fvCFD.H"
#include "polynomialMatrix.H"

namespace Foam
{
    class weightedMatrix : public polynomialMatrix
    {
    public:
        weightedMatrix(const scalarRectangularMatrix B);

        void applyStencilPointWeights(const scalarList& weights);
        void multiplyConstantAndLinearWeights(const scalar weight);
        void multiplyUpwindWeight(const scalar weight);
        void multiplyDownwindWeight(const scalar weight);
        scalarRectangularMatrix pseudoInverse() const;
    private:
        void multiplyStencilPointWeight(const label index, const scalar weight);
    };
}

#endif

