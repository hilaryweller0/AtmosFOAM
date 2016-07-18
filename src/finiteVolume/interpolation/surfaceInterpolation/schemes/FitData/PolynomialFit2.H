#ifndef PolynomialFit2_H
#define PolynomialFit2_H

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
    class PolynomialFit2
    {
    public:
        PolynomialFit2
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
        label numberOfSetBits(uint32_t i);

        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "PolynomialFit2.C"
#endif

#endif

