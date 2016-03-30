#ifndef fitCoefficients_H
#define fitCoefficients_H

#include "fvCFD.H"

namespace Foam
{
    class fitCoefficients
    {
    public:
        fitCoefficients
        (
            scalarList& coefficients,
            const List<point>& stencil,
            const bool linearCorrection,
            const scalar wLin
        );

        void applyCorrection(const bool goodFit);

    private:
        scalarList& coefficients;
        const bool linearCorrection;
        const scalar wLin;
    };
}

#endif

