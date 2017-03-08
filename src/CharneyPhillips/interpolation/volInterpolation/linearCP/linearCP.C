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
#include "linearCP.H"
#include "SVD.H"

template<class Type>
Foam::tmp<
    Foam::GeometricField
    <
        Type,
        Foam::fvPatchField,
        Foam::volMesh
    >
> Foam::linearCP<Type>::interpolate
(
    const Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>& s
) const
{
    return stencilDescription.weightedSum(s, coeffs);
}

template<class Type>
void Foam::linearCP<Type>::initCoeffs
(
    const fvMesh& mesh
)
{
    List<List<point> > stencilPoints(mesh.nCells());

    stencilDescription.collectData
    (
        mesh.Cf(),
        stencilPoints
    );

    forAll(stencilPoints, stencilForCellI)
    {
        const List<point>& stencil = stencilPoints[stencilForCellI];

        calculateInterpolationCoeffs
        (
            stencilDescription.stencil()[stencilForCellI],
            stencil,
            coeffs[stencilForCellI]
        );
    }
}

template<class Type>
void Foam::linearCP<Type>::calculateInterpolationCoeffs
(
    const labelList& stencilCellIndices,
    const List<point>& stencil,
    scalarList& coeffs
)
{
    // TODO
}
