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

        label size();
        scalar& operator[](int i);
        scalar operator[](int i) const;
        void applyCorrection(const bool goodFit);
        bool stable() const;

    private:
        scalarList& coefficients;
        const bool linearCorrection;
        const scalar wLin;
    };
}

#endif

