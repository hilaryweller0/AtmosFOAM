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

InNamespace
    Foam::fvc

Description
    Map from snGrad to gradient on edges

SourceFiles
    faceToPointReconstruct.C

\*---------------------------------------------------------------------------*/

#ifndef faceToPointReconstruct_H
#define faceToPointReconstruct_H

#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type, pointPatchField, pointMesh
        >
    > faceToPointReconstruct
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type, pointPatchField, pointMesh
        >
    > faceToPointReconstruct
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
    );
    
    // Apply coupled boundary conditions to point fields
    template<class Type>
    void applyCoupledBCs
    (
        GeometricField<Type, pointPatchField, pointMesh>& pf,
        const fvMesh& mesh
    );
    
    // Helper functions for applying coupled BCs to point fields
    // (taked from volPointInterpolate.C)
    template<class Type>
    void addSeparated(GeometricField<Type, pointPatchField, pointMesh>& pf);
    
    template<class Type>
    void pushUntransformedData
    (
        GeometricField<Type, pointPatchField, pointMesh>& pf,
        const fvMesh& mesh
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faceToPointReconstruct.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
