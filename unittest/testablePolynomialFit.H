#ifndef testablePolynomialFit_H
#define testablePolynomialFit_H

#include "fvCFD.H"
#include "fitCoefficients.H"
#include "fitResult.H"
#include "fitWeights.H"
#include "localStencil.H"

autoPtr<fitResult> fitPolynomial
(
        fitCoefficients& coefficients,
        fitWeights& weights,
        const localStencil& stencil
);

#endif
