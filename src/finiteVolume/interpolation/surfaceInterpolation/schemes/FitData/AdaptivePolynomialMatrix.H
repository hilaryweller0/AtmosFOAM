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
                const direction dimensions);

        scalarRectangularMatrix matrix();
    
    private:
        const List<point>& stencil;
        const direction dimensions;
        const label maxTerms;
    };
}

#ifdef NoRepository
#   include "AdaptivePolynomialMatrix.C"
#endif

#endif

