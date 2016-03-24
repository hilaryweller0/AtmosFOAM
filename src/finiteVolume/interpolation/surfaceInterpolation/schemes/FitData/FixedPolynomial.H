#ifndef FixedPolynomial_H
#define FixedPolynomial_H

#include "AbstractPolynomial.H"
#include "weightedMatrix.H"

namespace Foam
{
    template<class Polynomial>
    class FixedPolynomial : public AbstractPolynomial<Polynomial>
    {
    public:
        FixedPolynomial(
                const List<point>& stencil,
                const direction dimensions);

        scalarRectangularMatrix matrix() const;

    private:
        const List<point>& stencil;
        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "FixedPolynomial.C"
#endif

#endif

