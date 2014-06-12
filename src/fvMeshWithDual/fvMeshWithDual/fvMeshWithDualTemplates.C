/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "fvMeshWithDual.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::dualMap
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    if (sf.mesh() != *this)
    {
        FatalErrorIn("fvMeshWithDual::dualMap") << "cannot map " << sf.name()
            << " because it is on the wrong mesh" << exit(FatalError);
    }

//    Info << "Creating " << sf.name() << " on the dual mesh" << endl;
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(sf.name(), sf.instance(), dualMesh()),
            dualMesh(),
            dimensioned<Type>(sf.name(), sf.dimensions(), pTraits<Type>::zero),
            sf.boundaryField().types()
        )
    );
//    Info << "done" << endl;
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfd = tsfd();
    
    forAll(sf, faceI)
    {
        label faced = dualFaceMap()[faceI];
        if (faced != -1)
        {
            sfd[faced] = sf[faceI];
        }
    }
    forAll(sf.boundaryField(), patchI)
    {
        label faceI = sf.mesh().boundaryMesh()[patchI].start()-1;
        forAll(sf.boundaryField()[patchI], patchFace)
        {
            faceI++;
            label faced = dualFaceMap()[faceI];
            if (faced != -1)
            {
                if (faced < sfd.size())
                {
                    sfd[faced] = sf.boundaryField()[patchI][patchFace];
                }
                else
                {
                    const polyBoundaryMesh& bDualMesh=dualMesh().boundaryMesh();
                    const label dPatchI = bDualMesh.whichPatch(faced);
                    const label patchFaced = faced - bDualMesh[dPatchI].start();
                    sfd.boundaryField()[dPatchI][patchFaced]
                        = sf.boundaryField()[patchI][patchFace];
                }
            }
        }
    }
    
    return tsfd;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::dualMap
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        dualMap(tsf())
    );
    tsf.clear();
    return tsfd;

}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::signedDualMap
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(sf.name(), sf.instance(), dualMesh()),
            dualMesh(),
            dimensioned<Type>(sf.name(), sf.dimensions(), pTraits<Type>::zero),
            sf.boundaryField().types()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfd = tsfd();
    
    forAll(sf, faceI)
    {
        label faced = dualFaceMap()[faceI];
        if (faced != -1)
        {
            sfd[faced] = signMap(faceI)*sf[faceI];
        }
    }
    forAll(sf.boundaryField(), patchI)
    {
        label faceI = sf.mesh().boundaryMesh()[patchI].start()-1;
        forAll(sf.boundaryField()[patchI], patchFace)
        {
            faceI++;
            label faced = dualFaceMap()[faceI];
            if (faced != -1)
            {
                if (faced < sfd.size())
                {
                    sfd[faced] = signMap(faceI)
                                *sf.boundaryField()[patchI][patchFace];
                }
                else
                {
                    const polyBoundaryMesh& bDualMesh
                        =dualMesh().boundaryMesh();
                    const label dPatchI = bDualMesh.whichPatch(faced);
                    const label patchFaced = faced
                                           - bDualMesh[dPatchI].start();

                    sfd.boundaryField()[dPatchI][patchFaced]
                         = signMap(faceI)
                          *sf.boundaryField()[patchI][patchFace];
                }
            }
        }
    }
    
    return tsfd;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::signedDualMap
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        signedDualMap(tsf())
    );
    tsf.clear();
    return tsfd;

}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::dualFluxMap
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(sf.name(), sf.instance(), dualMesh()),
            dualMesh(),
            dimensioned<Type>(sf.name(), sf.dimensions(), pTraits<Type>::zero),
            sf.boundaryField().types()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfd = tsfd();
    
    forAll(sf, faceI)
    {
        label faced = dualFaceMap()[faceI];
        if (faced != -1)
        {
            sfd[faced] = signMap(faceI)*sf[faceI]*dualMesh().magSf()[faced]
                        /magSf()[faceI];
        }
    }
    forAll(sf.boundaryField(), patchI)
    {
        label faceI = sf.mesh().boundaryMesh()[patchI].start()-1;
        forAll(sf.boundaryField()[patchI], patchFace)
        {
            faceI++;
            label faced = dualFaceMap()[faceI];
            if (faced != -1)
            {
                if (faced < sfd.size())
                {
                    sfd[faced] = signMap(faceI)
                                *sf.boundaryField()[patchI][patchFace]
                                *dualMesh().magSf()[faced]
                                /magSf().boundaryField()[patchI][patchFace];
                }
                else
                {
                    const polyBoundaryMesh& bDualMesh
                        =dualMesh().boundaryMesh();
                    const label dPatchI = bDualMesh.whichPatch(faced);
                    const label patchFaced = faced
                                           - bDualMesh[dPatchI].start();

                    sfd.boundaryField()[dPatchI][patchFaced]
                         = signMap(faceI)
                          *sf.boundaryField()[patchI][patchFace]
                          *dualMesh().magSf().boundaryField()[dPatchI][patchFaced]
                          /magSf().boundaryField()[patchI][patchFace];
                }
            }
        }
    }
    
    return tsfd;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvMeshWithDual::dualFluxMap
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfd
    (
        dualFluxMap(tsf())
    );
    tsf.clear();
    return tsfd;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
