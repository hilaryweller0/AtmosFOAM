#ifndef PolynomialFit_H
#define PolynomialFit_H

#include "scalarList.H"

namespace Test {
    class PolynomialFit {
        public:
            PolynomialFit();
            const Foam::scalarList coefficients();

        private:
            Foam::scalarList coefficients_;
    };
}

#endif


