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

Class
    Foam::fvGhostMesh

Description
    Read in a mesh with ghost cells in order to handle periodic boundary
    conditions and provide functions for mapping between to the ghost cells

SourceFiles
    fvGhostMeshI.C
    fvGhostMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvGhostMesh_H
#define fvGhostMesh_H

#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fvGhostMesh Declaration
\*---------------------------------------------------------------------------*/

class fvGhostMesh
:
    public fvMesh
{
    // Private data

        //- Const reference to the original mesh
        const fvMesh& mesh_;
        
        //- The cell addressing from ghost cells onto the original cells
        IOList<label> ghostToMeshCells_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvGhostMesh(const fvGhostMesh&);

        //- Disallow default bitwise assignment
        void operator=(const fvGhostMesh&);


public:

    // Constructors

        //- Construct from components
        fvGhostMesh(const IOobject& io, const fvMesh& mesh__);


    //- Destructor
    ~fvGhostMesh();


    // Member Functions

        //- Map from a volField to the ghost volField
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > mapToGhost
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        //- Map from a volField to the ghost volField
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > mapToGhost
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >& vf
        )
        {
            return mapToGhost(vf());
            vf.clear();
        }

        //- Map from the ghost volField to the volField
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > mapFromGhost
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        //- Map from the ghost volField to the volField
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > mapFromGhost
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >& vf
        )
        {
            return mapFromGhost(vf());
            vf.clear();
        }

        //- Map from a surfaceField to the ghost surfaceField
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > mapToGhost
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );
        
        //- Map from a surfaceField to the ghost surfaceField
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > mapToGhost
        (
            const tmp<GeometricField<Type, fvPatchField, surfaceMesh> >& sf
        )
        {
            return mapToGhost(sf());
            sf.clear();
        }

        //- Map from the ghost surfaceField to the surfaceField
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > mapFromGhost
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );
        
        //- Map from the ghost surfaceField to the surfaceField
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > mapFromGhost
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& sf
        )
        {
            return mapFromGhost(sf());
            sf.clear();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvGhostMeshTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
