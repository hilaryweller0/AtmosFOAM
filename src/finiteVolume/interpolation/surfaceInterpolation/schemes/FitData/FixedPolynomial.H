#ifndef FixedPolynomial_H
#define FixedPolynomial_H

#include "AbstractPolynomial.H"

namespace Foam
{
    template<class Polynomial>
    class FixedPolynomial : public AbstractPolynomial<Polynomial>
    {
    public:
        FixedPolynomial(
                const List<point>& stencil,
                const direction dimensions);

        autoPtr<scalarRectangularMatrix> matrix() const;

    private:
        const List<point>& stencil;
        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "FixedPolynomial.C"
#endif

#endif

