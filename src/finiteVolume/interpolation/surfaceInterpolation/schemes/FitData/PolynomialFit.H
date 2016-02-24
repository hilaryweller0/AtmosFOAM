#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"

namespace Foam
{
    template<class Polynomial>
    class PolynomialFit
    {
    public:
        PolynomialFit
        (
            const bool linearCorrection,
            const scalar linearLimitFactor,
            const scalar centralWeight,
            const label dimension,
            const label terms
        );
        
        void fit
        (
            scalarList& coeffsi,
            scalarList& wts,
            const List<point>& stencilPoints,
            const scalar wLin,
            const point p0,
            const bool pureUpwind,
            const vector idir,
            const vector jdir,
            const vector kdir,
            const label faceI // only used for logging
        );
    
    private:
        const bool linearCorrection_;
        const scalar linearLimitFactor_;
        const scalar centralWeight_;
        const label dim_;
        const label minSize_;
    };
}

#ifdef NoRepository
#   include "PolynomialFit.C"
#endif

#endif
