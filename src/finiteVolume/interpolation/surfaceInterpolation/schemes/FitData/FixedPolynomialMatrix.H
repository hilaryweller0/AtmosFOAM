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
//    private:
        scalarRectangularMatrix B;
    };
}

#ifdef NoRepository
#   include "FixedPolynomialMatrix.C"
#endif

#endif

