/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef intersectEdge_H
#define intersectEdge_H

#include "PrimitivePatch.H"
#include "plane.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- The intersection on the sphere between edges ea and eb with radius on edge ea
inline point intersectEdge
(
    const primitivePatch& ePatch,
    const label ei,
    const point& a,
    const point& b
)
{
    const edge& e = ePatch.edges()[ei];
    const point& p = ePatch.localPoints()[e[0]];
    const point& q = ePatch.localPoints()[e[1]];
    const plane edgePlane(vector::zero, p, q);
    const plane abPlane(vector::zero, a, b);
    const plane::ray r = edgePlane.planeIntersect(abPlane);
    point pi = r.dir()/mag(r.dir())*sign(r.dir() & a);

    return pi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
