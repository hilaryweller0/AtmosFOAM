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
    const localStencil& stencil
)
{
    Polynomial polynomial(stencil, dimensions);
    autoPtr<scalarRectangularMatrix> B = polynomial.matrix();

    stabiliser stabiliser;
    const label columns = stabiliser.stabilise(B(), weights, coefficients, stencil);
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
