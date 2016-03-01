#ifndef AdaptivePolynomialMatrix_H
#define AdaptivePolynomialMatrix_H

namespace Foam
{
    template<class Polynomial>
    class AdaptivePolynomialMatrix
    {
    public:
        AdaptivePolynomialMatrix(
                const List<point>& stencil,
                const direction dimensions,
                const scalar tolerance=1e-4);

        scalarRectangularMatrix matrix();
    
    private:
        bool containsEntry(const List<label> fittableTerms, label term);
        bool fullRank(const scalarRectangularMatrix& B);

        const List<point>& stencil;
        const direction dimensions;
        const label maxTerms;
        const scalar tolerance;
    };
}

#ifdef NoRepository
#   include "AdaptivePolynomialMatrix.C"
#endif

#endif

