#ifndef TestPolynomialFit_H
#define TestPolynomialFit_H

#include "fvCFD.H"
#include "localStencil.H"

namespace Test {
    class PolynomialFit {
        public:
            PolynomialFit
            (
                const Foam::List<point>& stencilPoints,
                const Foam::label faceI
            );

            PolynomialFit
            (
                const Foam::localStencil stencil,
                bool linearCorrection=false
            );

            const Foam::scalarList coefficients();

        private:
            Foam::scalarList coefficients_;
    };
}

#endif


