#ifndef nullAdjuster_H
#define nullAdjuster_H

#include "fitCoefficients.H"
#include "weightedMatrix.H"

namespace Foam
{
    class nullAdjuster
    {
    public:
        nullAdjuster
        (
            const weightedMatrix& matrix, 
            fitCoefficients& coefficients,
            const fitWeights& weights
        );

        bool adjustWeights();

    private:
        const weightedMatrix& matrix;
        fitCoefficients& coefficients;
        const fitWeights& weights;
    };
}

#endif


