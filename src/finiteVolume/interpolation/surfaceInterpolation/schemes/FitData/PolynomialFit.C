#include "PolynomialFit.H"
#include "fitCoefficients.H"
#include "stabiliser.H"

template<class Polynomial>
Foam::PolynomialFit<Polynomial>::PolynomialFit
(
    const direction dimensions
)
:
    dimensions(dimensions)
{}

template<class Polynomial>
autoPtr<fitResult> Foam::PolynomialFit<Polynomial>::fit
(
    fitCoefficients& coefficients,
    fitWeights& weights,
    const localStencil& stencil,
    label faceI,
    bool owner
)
{
    Polynomial polynomial(stencil, dimensions);
    autoPtr<scalarRectangularMatrix> B = polynomial.matrix();

    stabiliser stabiliser;
    const label columns = stabiliser.stabilise(B(), weights, coefficients, faceI, owner);
    bool goodFit = (columns > 0);

    coefficients.applyCorrection(goodFit);

    return autoPtr<fitResult>(new fitResult(
            stencil,
            coefficients,
            weights,
            goodFit,
            columns
    ));
}
