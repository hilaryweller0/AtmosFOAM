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
            weightedMatrix& matrix, 
            scalarList& coefficients,
            scalarList& wts
        );

        bool adjustWeights();

    private:
        weightedMatrix& matrix;
        scalarList& coefficients;
        scalarList& wts;
    };
}

#endif


