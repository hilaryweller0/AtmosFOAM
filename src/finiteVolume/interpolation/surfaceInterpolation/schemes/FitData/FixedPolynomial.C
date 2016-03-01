#include "FixedPolynomial.H"

template<class Polynomial>
Foam::FixedPolynomial<Polynomial>::FixedPolynomial(
        const List<point>& stencil,
        const direction dimensions)
:
stencil(stencil),
dimensions(dimensions)
{}

template<class Polynomial>
scalarRectangularMatrix Foam::FixedPolynomial<Polynomial>::matrix() const
{
    scalarRectangularMatrix B(stencil.size(), Polynomial::nTerms(dimensions), scalar(0));

    forAll(stencil, i)
    {
        Polynomial::addCoeffs(B[i], stencil[i], 1, dimensions);
    }

    return B;
}

