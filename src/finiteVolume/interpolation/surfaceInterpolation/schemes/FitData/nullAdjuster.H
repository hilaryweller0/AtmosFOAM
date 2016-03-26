#ifndef nullAdjuster_H
#define nullAdjuster_H

#include "adjuster.H"
#include "polynomialMatrix.H"

namespace Foam
{
    class nullAdjuster : public adjuster
    {
    public:
        nullAdjuster
        (
            polynomialMatrix& matrix, 
            scalarList& coefficients,
            scalarList& wts
        );

        bool adjustWeights();

    private:
        polynomialMatrix& matrix;
        scalarList& coefficients;
        scalarList& wts;
    };
}

#endif


