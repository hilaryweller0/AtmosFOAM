#ifndef nullAdjuster_H
#define nullAdjuster_H

#include "adjuster.H"
#include "weightedMatrix.H"

namespace Foam
{
    class nullAdjuster : public adjuster
    {
    public:
        nullAdjuster
        (
            const weightedMatrix& matrix, 
            scalarList& coefficients,
            const fitWeights& weights
        );

        bool adjustWeights();

    private:
        const weightedMatrix& matrix;
        scalarList& coefficients;
        const fitWeights& weights;
    };
}

#endif


