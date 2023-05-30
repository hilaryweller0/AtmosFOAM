/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::weightedReconstructData

Description
    Data for weighted reconstructing cell vectors from face fluxes with
    specified weight for boundaries

SourceFiles
    weightedReconstructData.C

\*---------------------------------------------------------------------------*/

#ifndef weightedReconstructData_H
#define weightedReconstructData_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class weightedReconstructData Declaration
\*---------------------------------------------------------------------------*/

class weightedReconstructData
:
    public MeshObject<fvMesh, MoveableMeshObject, weightedReconstructData>
{
    // Private data

        //- The boundary weight
        const scalar boundaryWeight_;

        //- The weights for each face
        surfaceVectorField faceWeights_;

        //- The inverted tensor to convert from face flux to cell vectors
        volTensorField invTensor_;

    // Private Member Functions

        //- Re-calculate the data
        void calcData();

public:

    // Declare name of the class and its debug switch
    TypeName("weightedReconstructData");

    // Constructors

        //- Construct from components
        explicit weightedReconstructData
        (
            const fvMesh& mesh,
            const scalar boundaryWeight__
        );


    //- Destructor
    virtual ~weightedReconstructData();


    // Member functions

        //- Return the weights for each face (distance to weightPower)
        const surfaceVectorField& faceWeights() const { return faceWeights_;}

        //- Return he inverted tensor to convert from face flux to cell vectors
        const volTensorField& invTensor() const { return invTensor_;}

        //- Recalculate the data when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "weightedReconstructData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
