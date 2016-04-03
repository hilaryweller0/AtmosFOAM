#include "PolynomialFit.H"
#include "weightedMatrix.H"
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

    weightedMatrix matrix(B(), weights);

    stabiliser stabiliser;
    bool goodFit = stabiliser.stabilise(matrix, coefficients);

    coefficients.applyCorrection(goodFit);

    return autoPtr<fitResult>(new fitResult(
            stencil,
            coefficients,
            goodFit,
            B->m()
    ));
}
