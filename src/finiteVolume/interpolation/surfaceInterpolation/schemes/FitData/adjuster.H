#ifndef adjuster_H
#define adjuster_H

#include "polynomialMatrix.H"

namespace Foam
{
    class adjuster
    {
    public:
        adjuster
        (
            polynomialMatrix& matrix, 
            scalarList& coefficients,
            scalarList& wts,
            const scalar wLin,
            const bool pureUpwind,
            const bool linearCorrection,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );

        bool isGoodFit();

        void increaseWeights(bool firstIteration);

    private:
        bool eitherUpwindOrDownwindHasMaximumMagnitude();
        bool upwindCoefficientLargerThanSumOfOtherPositiveCoefficients();

        polynomialMatrix& matrix;
        scalarList& coefficients;
        scalarList& wts;
        const scalar wLin;
        const bool pureUpwind;
        const bool linearCorrection;
        const scalar linearLimitFactor;
        const scalar centralWeight;
    };
}

#endif

