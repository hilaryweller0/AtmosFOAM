#ifndef stabiliser_H
#define stabiliser_H

#include "fitCoefficients.H"
#include "fitWeights.H"

namespace Foam
{
    class stabiliser
    {
    public:
        label stabilise
        (
            const scalarRectangularMatrix& B,
            fitWeights& weights,
            fitCoefficients& coefficients
        ) const;
    };
}

#endif
