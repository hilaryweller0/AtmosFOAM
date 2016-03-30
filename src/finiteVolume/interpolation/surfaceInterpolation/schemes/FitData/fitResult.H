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
                const List<point>& stencilPoints,
                const fitCoefficients& coefficients,
                const bool good,
                const label polynomialTerms
        );

        const List<point>& stencilPoints;
        const fitCoefficients& coefficients;
        const bool good;
        const label polynomialTerms;
    };
}

#endif
