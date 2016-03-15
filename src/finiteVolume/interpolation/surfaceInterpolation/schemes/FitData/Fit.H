#ifndef Fit_H
#define Fit_H

#include "fvCFD.H"

namespace Foam
{
    class Fit
    {
    public:
        Fit(
                const List<point>& stencilPoints,
                const scalarList& coeffs,
                const bool good
        );

        const List<point>& stencilPoints;
        const scalarList& coeffs;
        const bool good;
    };
}

#endif
