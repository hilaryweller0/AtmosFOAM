/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fvcLocalMinMax.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void localMin
(
    Field<scalar>& vfMin,
    const surfaceScalarField& ssf
)
{
    const fvMesh& mesh = ssf.mesh();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vfMin[owner[facei]] = min(vfMin[owner[facei]], ssf[facei]);
        vfMin[neighbour[facei]] = min(vfMin[neighbour[facei]], ssf[facei]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<scalar>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vfMin[pFaceCells[facei]] = min(vfMin[pFaceCells[facei]], pssf[facei]);
        }
    }
}


tmp<volScalarField>
localMin
(
    const surfaceScalarField& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<volScalarField> tvf
    (
        new volScalarField
        (
            IOobject
            (
                "localMin("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>
            (
                "GREAT",
                ssf.dimensions(),
                GREAT
            ),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
    volScalarField& vf = tvf.ref();

    localMin(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


tmp<volScalarField>
localMin
(
    const tmp<surfaceScalarField>& tssf
)
{
    tmp<volScalarField> tvf
    (
        localMin(tssf())
    );
    tssf.clear();
    return tvf;
}


void localMax
(
    Field<scalar>& vfMax,
    const surfaceScalarField& ssf
)
{
    const fvMesh& mesh = ssf.mesh();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vfMax[owner[facei]] = max(vfMax[owner[facei]], ssf[facei]);
        vfMax[neighbour[facei]] = max(vfMax[neighbour[facei]], ssf[facei]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<scalar>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vfMax[pFaceCells[facei]] = max(vfMax[pFaceCells[facei]], pssf[facei]);
        }
    }
}


tmp<volScalarField>
localMax
(
    const surfaceScalarField& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<volScalarField> tvf
    (
        new volScalarField
        (
            IOobject
            (
                "localMax("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>
            (
                "-GREAT",
                ssf.dimensions(),
                -GREAT
            ),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
    volScalarField& vf = tvf.ref();

    localMax(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


tmp<volScalarField>
localMax
(
    const tmp<surfaceScalarField>& tssf
)
{
    tmp<volScalarField> tvf
    (
        localMax(tssf())
    );
    tssf.clear();
    return tvf;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
