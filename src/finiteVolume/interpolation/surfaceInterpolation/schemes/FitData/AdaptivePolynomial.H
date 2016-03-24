#ifndef AdaptivePolynomial_H
#define AdaptivePolynomial_H

#include "AbstractPolynomial.H"

namespace Foam
{
    /**
     * Calculate the largest polynomial matrix that is full rank.
     * If there are more terms than stencil points, or some stencil
     * points are very closely spaced, then the lowest order terms
     * will be included and the higher order terms dropped.
     */
    template<class Polynomial>
    class AdaptivePolynomial : public AbstractPolynomial<Polynomial>
    {
    public:
        AdaptivePolynomial(
                const List<point>& stencil,
                const direction dimensions,
                const scalar tolerance=1e-4);

        autoPtr<scalarRectangularMatrix> matrix() const;
    
    private:
        bool containsEntry(const List<label> fittableTerms, label term) const;
        bool fullRank(const scalarRectangularMatrix& B) const;

        const List<point>& stencil;
        const direction dimensions;
        const label maxTerms;
        const scalar tolerance;
    };
}

#ifdef NoRepository
#   include "AdaptivePolynomial.C"
#endif

#endif

