#ifndef fitResult_H
#define fitResult_H

#include "fvCFD.H"
#include "fitCoefficients.H"
#include "fitWeights.H"
#include "localStencil.H"

namespace Foam
{
    class fitResult
    {
    public:
        fitResult
        (
                const localStencil stencil,
                const fitCoefficients& coefficients,
                const fitWeights& weights,
                const bool good,
                const label polynomialTerms
        );

        const localStencil stencil;
        const fitCoefficients& coefficients;
        const fitWeights weights;
        const bool good;
        const label polynomialTerms;
    };
}

#endif
