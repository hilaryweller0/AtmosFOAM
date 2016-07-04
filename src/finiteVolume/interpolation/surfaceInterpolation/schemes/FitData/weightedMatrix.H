#ifndef weightedMatrix_H
#define weightedMatrix_H

#include <limits>
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
        template<typename _Matrix_Type_> _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon()) const;

        scalarRectangularMatrix B;
        const autoPtr<fitWeights> weights;
    };
}

#endif

