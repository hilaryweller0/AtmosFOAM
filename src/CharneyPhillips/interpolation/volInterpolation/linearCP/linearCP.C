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
    Foam::tmp<
        Foam::GeometricField
        <
            Type,
            Foam::fvPatchField,
            Foam::volMesh
        >
    > tInterpField = stencilDescription.weightedSum(s, coeffs);
    Foam::GeometricField
    <
        Type,
        Foam::fvPatchField,
        Foam::volMesh
    >& interpField = tInterpField.ref();

    // copy boundary values
    forAll(s.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& sourcePatch = s.boundaryField()[patchI];
        fvPatchField<Type>& targetPatch = interpField.boundaryFieldRef()[patchI];
        targetPatch = sourcePatch;
    }

    return tInterpField;
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
            mesh.C()[stencilForCellI],
            stencilDescription.stencil()[stencilForCellI],
            stencil,
            coeffs[stencilForCellI]
        );
    }
}

template<class Type>
void Foam::linearCP<Type>::calculateInterpolationCoeffs
(
    const point& origin,
    const labelList& stencilFaceIndices,
    const List<point>& stencil,
    scalarList& coeffs
)
{
    label nonVerticalFaces = 0;
    boolList include(stencilFaceIndices.size());
    forAll(stencilFaceIndices, i)
    {
        // non-vertical faces are included
        include[i] = mag(g.unitFaceNormal()[stencilFaceIndices[i]]) > 1e-12;
        if (include[i]) nonVerticalFaces++;
    }

    scalarRectangularMatrix B(nonVerticalFaces, 2);

    label n = 0;
    forAll(stencil, i)
    {
        if (include[i])
        {
            const point& p = stencil[i] - origin;

            B(n, 0) = 1;
            B(n, 1) = p.z();

            n++;
        }
    }

    const SVD svd(B);
    const scalarRectangularMatrix& Binv = svd.VSinvUt();

    n = 0;
    forAll(stencil, i)
    {
        if (include[i])
        {            
            coeffs.append(Binv(0, n));
            n++;
        }
        else
        {
            coeffs.append(0);
        }
    }
}
