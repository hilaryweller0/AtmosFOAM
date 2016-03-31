#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "Basis.H"
#include "fitResult.H"
#include "weightedMatrix.H"
#include "fitWeights.H"
#include "localStencil.H"

namespace Foam
{
    template<class Polynomial>
    class PolynomialFit
    {
    public:
        PolynomialFit
        (
            const direction dimensions
        );
        
        autoPtr<fitResult> fit
        (
            fitCoefficients& coefficients,
            fitWeights& weights,
            const localStencil& stencil
        );
    
    private:
        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "PolynomialFit.C"
#endif

#endif
