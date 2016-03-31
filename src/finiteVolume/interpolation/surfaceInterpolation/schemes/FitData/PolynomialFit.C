#include "PolynomialFit.H"
#include "weightedMatrix.H"
#include "nullAdjuster.H"
#include "fitCoefficients.H"

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
    weightedMatrix matrix(B());

    matrix.apply(weights);

    nullAdjuster adjuster
    (
        matrix, 
        coefficients,
        weights
    );

    bool goodFit = adjuster.adjustWeights();
    coefficients.applyCorrection(goodFit);

    return autoPtr<fitResult>(new fitResult(
            coefficients,
            goodFit,
            B->m()
    ));
}
