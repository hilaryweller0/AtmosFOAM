#ifndef weightedMatrix_H
#define weightedMatrix_H

#include "fvCFD.H"
#include "fitWeights.H"

namespace Foam
{
    class weightedMatrix
    {
    public:
        weightedMatrix(scalarRectangularMatrix& B);

        void apply(const fitWeights& weights);
        scalarRectangularMatrix pseudoInverse() const;

    private:
        void applyStencilPointWeights(const fitWeights& weights);

        scalarRectangularMatrix& B;
    };
}

#endif

