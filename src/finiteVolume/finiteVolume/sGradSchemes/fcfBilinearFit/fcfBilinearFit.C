/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fvc.H"
#include "fcfBilinearFit.H"

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvsPatchField,
        Foam::surfaceMesh
    >
> Foam::fv::fcfBilinearFit<Type>::operator()
(
    const Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>& s
) const
{
//    return fvc::interpolate(fvc::grad(s));
    return stencil.weightedSum(s, coeffs);
}

template<class Type>
void Foam::fv::fcfBilinearFit<Type>::initCoeffs
(
    const fvMesh& mesh
)
{
    List<List<point> > stencilPoints(mesh.nFaces());

    stencil.collectData
    (
        mesh.Cf(),
        stencilPoints
    );

    forAll(stencilPoints, stencilForFaceI)
    {
        // TODO: construct a matrix using the stencil points and get the pseudo-inverse
        forAll(stencilPoints[stencilForFaceI], faceIInStencil)
        {
            coeffs[stencilForFaceI].append(vector(1, 0, 0));
        }
    }
}
