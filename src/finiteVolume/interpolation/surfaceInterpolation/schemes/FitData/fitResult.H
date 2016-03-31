#ifndef fitResult_H
#define fitResult_H

#include "fvCFD.H"
#include "fitCoefficients.H"

namespace Foam
{
    class fitResult
    {
    public:
        fitResult(
                const fitCoefficients& coefficients,
                const bool good,
                const label polynomialTerms
        );

        const fitCoefficients& coefficients;
        const bool good;
        const label polynomialTerms;
    };
}

#endif
