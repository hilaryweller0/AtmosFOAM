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
#include "fvcReconstructCP.H"

template<class Type>
Foam::tmp<
    Foam::GeometricField
    <
        Type,
        Foam::fvPatchField,
        Foam::volMesh
    >
> Foam::fvcReconstructCP<Type>::interpolate
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
    > tInterpField = fvc::reconstruct(
            s
            * g.unitFaceNormal()
            * volInterpolationScheme<Type>::mesh().magSf()
        ) & g.unit();
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
