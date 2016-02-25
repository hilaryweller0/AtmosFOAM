#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"
#include "Basis.H"

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
            const direction dimensions,
            const label terms
        );
        
        void fit
        (
            scalarList& coeffsi,
            scalarList& wts,
            const List<point>& stencilPoints,
            const scalar wLin,
            const point& p0,
            const bool pureUpwind,
            const Basis& basis,
            const label faceI // only used for logging
        );
    
    private:
        vector toLocalCoordinates(
                const point& origin,
                const point& p,
                const Basis& basis);

        const bool linearCorrection_;
        const scalar linearLimitFactor_;
        const scalar centralWeight_;
        const direction dim_;
        const label minSize_;
    };
}

#ifdef NoRepository
#   include "PolynomialFit.C"
#endif

#endif
