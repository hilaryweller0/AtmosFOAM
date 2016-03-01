#ifndef FixedPolynomialMatrix_H
#define FixedPolynomialMatrix_H

#include "PolynomialMatrix.H"

namespace Foam
{
    template<class Polynomial>
    class FixedPolynomialMatrix : public PolynomialMatrix<Polynomial>
    {
    public:
        FixedPolynomialMatrix(
                const List<point>& stencil,
                const direction dimensions);

        void applyStencilPointWeights(const scalarList& weights);
        void multiplyConstantAndLinearWeights(const scalar weight);
        void multiplyUpwindWeight(const scalar weight);
        void multiplyDownwindWeight(const scalar weight);
        void multiplyStencilPointWeight(const label index, const scalar weight);
        scalarRectangularMatrix matrix() const;
    private:
        scalarRectangularMatrix B;
    };
}

#ifdef NoRepository
#   include "FixedPolynomialMatrix.C"
#endif

#endif

