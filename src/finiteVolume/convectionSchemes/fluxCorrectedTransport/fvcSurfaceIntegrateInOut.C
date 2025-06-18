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

#include "fvcSurfaceIntegrateInOut.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrateIn
(
    Field<Type>& ivf,
    const SurfaceField<Type>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ivf[owner[facei]] -= min(ssf[facei],0);
        ivf[neighbour[facei]] += max(ssf[facei], 0);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] -= min(pssf[facei],0);
        }
    }

    ivf /= mesh.Vsc();
}


template<class Type>
tmp<VolField<Type>>
surfaceIntegrateIn
(
    const SurfaceField<Type>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> tvf
    (
        new VolField<Type>
        (
            IOobject
            (
                "surfaceIntegrateIn("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimVolume,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );
    VolField<Type>& vf = tvf.ref();

    surfaceIntegrateIn(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<VolField<Type>>
surfaceIntegrateIn
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    tmp<VolField<Type>> tvf
    (
        surfaceIntegrateIn(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
void surfaceIntegrateOut
(
    Field<Type>& ivf,
    const SurfaceField<Type>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ivf[owner[facei]] += max(ssf[facei],0);
        ivf[neighbour[facei]] -= min(ssf[facei], 0);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += max(pssf[facei],0);
        }
    }

    ivf /= mesh.Vsc();
}


template<class Type>
tmp<VolField<Type>>
surfaceIntegrateOut
(
    const SurfaceField<Type>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<VolField<Type>> tvf
    (
        new VolField<Type>
        (
            IOobject
            (
                "surfaceIntegrateOut("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimVolume,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );
    VolField<Type>& vf = tvf.ref();

    surfaceIntegrateOut(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<VolField<Type>>
surfaceIntegrateOut
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    tmp<VolField<Type>> tvf
    (
        surfaceIntegrateOut(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
