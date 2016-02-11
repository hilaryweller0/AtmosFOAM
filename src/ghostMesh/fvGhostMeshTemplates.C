/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::fvGhostMesh::mapToGhost
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // Check that vf is on mesh_
    if (vf.mesh() != mesh_)
    {
        FatalErrorIn("fvGhostMesh::mapToGhost")
            << vf.name() << " is not on the original mesh "
            << exit(FatalError);
    }
    
    // Create the new vol field on the ghost mesh
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvfg
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject(vf.name(), mesh_.time().timeName(), *this),
            *this,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero),
            "zeroGradient"
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vfg = tvfg();
    
    // Transfer the original cells
    forAll(vf, cellI) {vfg[cellI] = vf[cellI];}
    
    // Transfers to ghost cells
    const label nCells0 = vf.mesh().nCells();
    forAll(ghostToMeshCells_, cellI)
    {
        vfg[nCells0+cellI] = vf[ghostToMeshCells_[cellI]];
    }
    
    return tvfg;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::fvGhostMesh::mapFromGhost
(
    const GeometricField<Type, fvPatchField, volMesh>& vfg
)
{
    // Create the new vol field on the original mesh
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject(vfg.name(), mesh_.time().timeName(), mesh_),
            mesh_,
            dimensioned<Type>(vfg.name(), vfg.dimensions(), pTraits<Type>::zero),
            "zeroGradient"
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    // Transfer the original cells
    forAll(vf, cellI) { vf[cellI] = vfg[cellI]; }
    
    return tvf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fvGhostMesh::mapToGhost
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    // Check that sf is on mesh_
    if (sf.mesh() != mesh_)
    {
        FatalErrorIn("fvGhostMesh::mapToGhost")
            << sf.name() << " is not on the original mesh "
            << exit(FatalError);
    }

    // Create the new surface field on the ghost mesh
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfg
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(sf.name()+"ghost", mesh_.time().timeName(), *this),
            *this,
            sf.dimensions(),
            sf.internalField(),
            sf.boundaryField()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfg = tsfg();
    
    return tsfg;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fvGhostMesh::mapFromGhost
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfg
)
{
    // Create the new surface field on the original mesh
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(sfg.name()+"FromGhost", mesh_.time().timeName(), mesh_),
            mesh_,
            sfg.dimensions(),
            sfg.internalField(),
            sfg.boundaryField()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();
    
    return tsf;
}

// ************************************************************************* //
