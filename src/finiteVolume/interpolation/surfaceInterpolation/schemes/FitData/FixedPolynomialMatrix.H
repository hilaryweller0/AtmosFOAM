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
        void multiplyConstantAndLinearWeights(const scalar weight);
        void multiplyUpwindWeight(const scalar weight);
        void multiplyDownwindWeight(const scalar weight);
        scalarRectangularMatrix pseudoInverse();
    private:
        scalarRectangularMatrix B;
        const direction dimensions;
    };
}

#ifdef NoRepository
#   include "FixedPolynomialMatrix.C"
#endif

#endif

