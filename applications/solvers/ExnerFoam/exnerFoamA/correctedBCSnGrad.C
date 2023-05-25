/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "correctedBCSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::correctedBCSnGrad<Type>::~correctedBCSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::correctedBCSnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    > grad
    (
        gradScheme<Type>::New
        (
            mesh,
            mesh.schemes().grad("grad(" + vf.name() + ')')
        )().grad(vf, "grad(" + vf.name() + ')')
    );

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf =
        linear<typename outerProduct<vector, Type>::type>(mesh).dotInterpolate
        (
            mesh.nonOrthCorrectionVectors(),
            grad
        );
    tssf.ref().rename("snGradCorr(" + vf.name() + ')');

    // Extrapolation correction at boundaries
    GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> unCorrGrad
        = snGradScheme<Type>::snGrad(vf, deltaCoeffs(vf));
    
    forAll(tssf.ref().boundaryFieldRef(), patchi)
    {
        fvsPatchField<Type>& bCorr = tssf.ref().boundaryFieldRef()[patchi];
        vectorField Sfhat = mesh.Sf().boundaryField()[patchi]
                            / mesh.magSf().boundaryField()[patchi];
        if(!bCorr.coupled())
        {
            bCorr = unCorrGrad.boundaryFieldRef()[patchi]
                  - (grad.boundaryField()[patchi].patchInternalField() & Sfhat);
        }
    }

    return tssf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::correctedBCSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "snGradCorr("+vf.name()+')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        ssf.replace
        (
            cmpt,
            correctedBCSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
