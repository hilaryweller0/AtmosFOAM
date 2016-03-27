#ifndef iterativeAdjuster_H
#define iterativeAdjuster_H

#include "adjuster.H"
#include "weightedMatrix.H"

namespace Foam
{
    class iterativeAdjuster : public adjuster
    {
    public:
        iterativeAdjuster
        (
            weightedMatrix& matrix, 
            scalarList& coefficients,
            fitWeights& weights,
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

        weightedMatrix& matrix;
        scalarList& coefficients;
        fitWeights& weights;
        const scalar wLin;
        const bool pureUpwind;
        const bool linearCorrection;
        const scalar linearLimitFactor;
        const scalar centralWeight;
    };
}

#endif

