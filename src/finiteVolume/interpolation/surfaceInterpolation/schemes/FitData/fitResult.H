#ifndef fitResult_H
#define fitResult_H

#include "fvCFD.H"

namespace Foam
{
    class fitResult
    {
    public:
        fitResult(
                const List<point>& stencilPoints,
                const scalarList& coeffs,
                const bool good,
                const label polynomialTerms
        );

        const List<point>& stencilPoints;
        const scalarList& coeffs;
        const bool good;
        const label polynomialTerms;
    };
}

#endif
