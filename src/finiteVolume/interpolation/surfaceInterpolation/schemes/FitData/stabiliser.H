#ifndef stabiliser_H
#define stabiliser_H

#include "fitCoefficients.H"
#include "fitWeights.H"
#include "weightedMatrix.H"

namespace Foam
{
    class stabiliser
    {
    public:
        label stabilise
        (
            const scalarRectangularMatrix& B,
            fitWeights& weights,
            fitCoefficients& coefficients,
            label faceI,
            bool owner
        ) const;

    private:
        autoPtr<weightedMatrix> findStabilisableMatrix
        (
            const scalarRectangularMatrix& B,
            const fitWeights& weights,
            const fitCoefficients& coefficients,
            label faceI,
            bool owner
        ) const;
    };
}

#endif
