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
        weightedMatrix
        (
            const scalarRectangularMatrix& matrix
        );

        weightedMatrix
        (
            const scalarRectangularMatrix& matrix,
            const fitWeights& weights,
            bool applyWeights = true
        );

        weightedMatrix
        (
            const weightedMatrix& matrix,
            const fitWeights& weights
        );

        void populate(fitCoefficients& coefficients) const;

        autoPtr<weightedMatrix> subset(const labelList columns) const;

        label columns() const;

    private:
        void applyStencilPointWeights(const fitWeights& weights);

        scalarRectangularMatrix B;
        const autoPtr<fitWeights> weights;
    };
}

#endif

