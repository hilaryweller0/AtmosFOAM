#ifndef iterativeAdjuster_H
#define iterativeAdjuster_H

#include "adjuster.H"
#include "polynomialMatrix.H"

namespace Foam
{
    class iterativeAdjuster : public adjuster
    {
    public:
        iterativeAdjuster
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

        bool adjustWeights();

    private:
        bool isGoodFit();
        bool eitherUpwindOrDownwindHasMaximumMagnitude();
        bool upwindCoefficientLargerThanSumOfOtherPositiveCoefficients();
        void increaseWeights(bool firstIteration);

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

