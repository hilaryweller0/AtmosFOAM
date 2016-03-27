#ifndef weightedMatrix_H
#define weightedMatrix_H

#include "fvCFD.H"
#include "fitWeights.H"

namespace Foam
{
    class weightedMatrix
    {
    public:
        weightedMatrix(scalarRectangularMatrix& B);

        void apply(const fitWeights& weights);
        void multiplyConstantAndLinearWeights(const scalar weight); 
        void multiplyUpwindWeight(const scalar weight);
        void multiplyDownwindWeight(const scalar weight);
        scalarRectangularMatrix pseudoInverse() const;

    private:
        void applyStencilPointWeights(const fitWeights& weights);
        void multiplyStencilPointWeight(const label index, const scalar weight);

        scalarRectangularMatrix& B;
    };
}

#endif

