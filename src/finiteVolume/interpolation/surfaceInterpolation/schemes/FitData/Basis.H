#ifndef Basis_H
#define Basis_H

#include "fvCFD.H"

namespace Foam
{
    class Basis
    {
        public:
            Basis(vector i, vector j, vector k);
            const vector i, j, k;
    };
}

#ifdef NoRepository
#   include "Basis.C"
#endif

#endif
