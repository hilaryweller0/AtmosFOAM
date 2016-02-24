#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "fvCFD.H"

namespace Test {
    class PolynomialFit {
        public:
            PolynomialFit(
                    const Foam::List<point>& stencilPoints,
                    const Foam::label faceI);
            const Foam::scalarList coefficients();

        private:
            Foam::scalarList coefficients_;
    };
}

#endif


