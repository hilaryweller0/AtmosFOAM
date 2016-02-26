#ifndef FixedPolynomialMatrix_H
#define FixedPolynomialMatrix_H

namespace Foam
{
    template<class Polynomial>
    class FixedPolynomialMatrix
    {
    public:
        FixedPolynomialMatrix(
                const List<point>& stencil,
                const direction dimensions);

        void setStencilPoint(const label index, const point& p);
        void applyStencilPointWeights(const scalarList& weights);
//    private:
        scalarRectangularMatrix B;
    private:
        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "FixedPolynomialMatrix.C"
#endif

#endif

