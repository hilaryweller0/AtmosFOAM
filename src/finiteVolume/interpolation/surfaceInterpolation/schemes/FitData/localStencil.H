#ifndef localStencil_H
#define localStencil_H

#include "fvCFD.H"
#include "Basis.H"

namespace Foam
{
    class localStencil
    {
    public:
        localStencil(const List<point> localPoints);

        localStencil(
            const List<point>& stencilPoints,
            const point& origin,
            const Basis& basis
        );

        void upwindUpwindIndices(List<label>& indices) const;
        label size() const;
        point operator[](int i) const;

    private:
        scalar scaleLocalCoordinates
        (
            const point& origin,
            const point& upwindPoint,
            const Basis& basis
        );

        point toLocalCoordinates
        (
            const point& origin,
            const point& p,
            const Basis& basis
        );

        List<point> localPoints;
    };
}

#endif
