#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"
#include "Basis.H"
#include "WeightedMatrix.H"

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
        
        void fit
        (
            scalarList& coeffsi,
            scalarList& wts,
            const List<point>& stencilPoints,
            const scalar wLin, // only used when linearCorrection=true
            const point& origin,
            bool pureUpwind,
            const Basis& basis,
            const label faceI // only used for logging
        );
    
    private:
        List<point> toLocalCoordinates(
                const List<point>& stencilPoints,
                const point& origin,
                const Basis& basis);

        point toLocalCoordinates(
                const point& origin,
                const point& p,
                const Basis& basis);

        scalar scaleLocalCoordinates(
                const point& origin,
                const point& upwindPoint,
                const Basis& basis);

        bool isGoodFit(const scalarList& coefficients, const scalar wLin);

        bool eitherUpwindOrDownwindHasMaximumMagnitude(
                const scalarList& coefficients);

        bool upwindCoefficientLargerThanSumOfOtherPositiveCoefficients(
                const scalarList& coefficients);

        void increaseWeights(
                WeightedMatrix& matrix,
                scalarList& wts,
                bool pureUpwind,
                bool firstIteration);

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
