#include "localStencil.H"

Foam::localStencil::localStencil(const List<point> localPoints)
:
    localPoints(localPoints)
{}

Foam::localStencil::localStencil(
    const List<point>& stencilPoints,
    const point& origin,
    const Basis& basis
)
:
localPoints(stencilPoints.size(), point(0, 0, 0))
{
    scalar scale = scaleLocalCoordinates(origin, stencilPoints[0], basis);
    forAll(stencilPoints, i)
    {
        localPoints[i] = toLocalCoordinates(origin, stencilPoints[i], basis) / scale;
    }
}

scalar Foam::localStencil::scaleLocalCoordinates
(
    const point& origin,
    const point& upwindPoint,
    const Basis& basis
)
{
    return cmptMax(cmptMag((toLocalCoordinates(origin, upwindPoint, basis))));
}

point Foam::localStencil::toLocalCoordinates
(
    const point& origin,
    const point& p,
    const Basis& basis
)
{
    point d;

    d.x() = (p - origin)&basis.i;
    d.y() = (p - origin)&basis.j;
    #ifndef SPHERICAL_GEOMETRY
    d.z() = (p - origin)&basis.k;
    #else
    d.z() = mag(p) - mag(origin);
    #endif

    return d;
}

point Foam::localStencil::operator[](int i) const
{
    return localPoints[i];
}

label Foam::localStencil::size() const
{
    return localPoints.size();
}
