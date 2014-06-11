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
    Foam::faceToFaceStencil

Description
    baseclass for extended face centred addressing. Contains per face a
    list of neighbouring faces in global addressing.

SourceFiles
    faceToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef faceToFaceStencil_H
#define faceToFaceStencil_H

#include "globalIndex.H"
#include "globalIndexAndTransform.H"
#include "boolList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class faceToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class faceToFaceStencil
{
    // Private data

        const polyMesh& mesh_;

        //- Global numbering for untransformed faces
        const globalIndex globalNumbering_;

        //- Global numbering for transformed faces
        const globalIndexAndTransform globalTransforms_;


protected:

    // Protected data

        //- Untransformed faces
        labelListList elements_;

        //- Transformed faces
        List<labelPairList> transformedElements_;


public:

    // Constructors

        //- Construct from mesh
        explicit faceToFaceStencil(const polyMesh&);


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Global numbering for untransformed faces
        const globalIndex& globalNumbering() const
        {
            return globalNumbering_;
        }

        //- Global numbering for transformed faces
        const globalIndexAndTransform& globalTransforms() const
        {
            return globalTransforms_;
        }

        //- Untransformed faces
        const labelListList& elements() const
        {
            return elements_;
        }

        //- Untransformed faces
        labelListList& elements()
        {
            return elements_;
        }

        //- Transformed faces
        const List<labelPairList>& transformedElements() const
        {
            return transformedElements_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
