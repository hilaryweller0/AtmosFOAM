#include "FixedPolynomialMatrix.H"

template<class Polynomial>
Foam::FixedPolynomialMatrix<Polynomial>::FixedPolynomialMatrix(
        const List<point>& stencil,
        const direction dimensions)
:
    B(stencil.size(), Polynomial::nTerms(dimensions), scalar(0))
{
}
