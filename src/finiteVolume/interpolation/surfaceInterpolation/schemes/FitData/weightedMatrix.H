#ifndef weightedMatrix_H
#define weightedMatrix_H

#include "fvCFD.H"
#include "fitCoefficients.H"
#include "fitWeights.H"

namespace Foam
{
    class weightedMatrix
    {
    public:
        weightedMatrix(
            scalarRectangularMatrix& B,
            const fitWeights& weights
        );

        void populate(fitCoefficients& coefficients) const;

    private:
        void applyStencilPointWeights(const fitWeights& weights);

        scalarRectangularMatrix& B;
        const fitWeights& weights;
    };
}

#endif

