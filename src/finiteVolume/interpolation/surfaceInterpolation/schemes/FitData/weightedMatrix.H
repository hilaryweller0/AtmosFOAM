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
            const scalarRectangularMatrix& matrix,
            const fitWeights& weights,
            bool applyWeights = true
        );

        void populate(fitCoefficients& coefficients) const;

        autoPtr<weightedMatrix> truncateTo(const label columns) const;

        label columns() const;

    private:
        void applyStencilPointWeights(const fitWeights& weights);

        scalarRectangularMatrix B;
        const fitWeights& weights;
    };
}

#endif

