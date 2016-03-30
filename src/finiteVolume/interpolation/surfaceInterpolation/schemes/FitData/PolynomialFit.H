#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "Basis.H"
#include "Fit.H"
#include "weightedMatrix.H"
#include "fitWeights.H"

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
            const direction dimensions
        );
        
        autoPtr<Fit> fit
        (
            scalarList& coeffsi,
            fitWeights& weights,
            const List<point>& stencilPoints,
            const scalar wLin, // only used when linearCorrection=true
            const point& origin,
            const bool pureUpwind,
            const Basis& basis
        );
    
    private:
        List<point> toLocalCoordinates
        (
            const List<point>& stencilPoints,
            const point& origin,
            const Basis& basis
        );

        point toLocalCoordinates
        (
            const point& origin,
            const point& p,
            const Basis& basis
        );

        scalar scaleLocalCoordinates
        (
            const point& origin,
            const point& upwindPoint,
            const Basis& basis
        );

        const bool linearCorrection_;
        const scalar linearLimitFactor_;
        const scalar centralWeight_;
        const direction dim_;
    };
}

#ifdef NoRepository
#   include "PolynomialFit.C"
#endif

#endif
