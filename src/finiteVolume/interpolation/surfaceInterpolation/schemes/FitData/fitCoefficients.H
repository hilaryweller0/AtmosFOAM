#ifndef fitCoefficients_H
#define fitCoefficients_H

#include "fvCFD.H"
#include "localStencil.H"

namespace Foam
{
    class fitCoefficients
    {
    public:
        fitCoefficients(const fitCoefficients& c);

        fitCoefficients
        (
            const label stencilSize,
            const bool linearCorrection,
            const scalar wLin
        );

        label size() const;
        scalar& operator[](int i);
        scalar operator[](int i) const;
        friend Ostream& operator<< (Ostream& stream, const fitCoefficients& c);
        void copyFrom(const scalarList& source);
        void copyFrom(const fitCoefficients& source);
        void copyInto(scalarList& target);
        void applyCorrection(const bool goodFit);
        scalar sum() const;

    private:
        scalarList coefficients;
        const bool linearCorrection;
        const scalar wLin;
    };
}

#endif

