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
    Foam::extendedCentredFaceToFaceStencil

Description
    Stencil to get access to faces on owner and neighbour cells.

SourceFiles
    extendedCentredFaceToFaceStencil.C
    extendedCentredFaceToFaceStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredFaceToFaceStencil_H
#define extendedCentredFaceToFaceStencil_H

#include "labelList.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceToFaceStencil;
class mapDistribute;

/*---------------------------------------------------------------------------*\
              Class extendedCentredFaceToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredFaceToFaceStencil
{
    // Private data

        //- Swap map for getting neigbouring data
        autoPtr<mapDistribute> mapPtr_;

        //- Per face the untransformed elements
        labelListList elements_;

        //- Per face the transformed elements
        labelListList transformedElements_;


    // Private Member Functions

        //- Use map to get the (non-position) data into stencil order
        template<class T, class TransformOp>
        void collectData
        (
            const GeometricField<T, fvsPatchField, surfaceMesh>& fld,
            List<List<T> >& stencilFld,
            const TransformOp& top
        ) const;

        //- Disallow default bitwise copy construct
        extendedCentredFaceToFaceStencil
        (
            const extendedCentredFaceToFaceStencil&
        );

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredFaceToFaceStencil&);


public:

    // Constructors

        //- Construct from uncompacted face stencil
        explicit extendedCentredFaceToFaceStencil(const faceToFaceStencil&);


    // Member Functions

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            return mapPtr_();
        }

        //- Return reference to the untransformed elements in the stencil
        const labelListList& elements() const
        {
            return elements_;
        }

        //- Return reference to the transformed elements in the stencil
        const labelListList& transformedElements() const
        {
            return transformedElements_;
        }

        //- After removing elements from the stencil adapt the schedule (map).
        void compact();

        //- Use map to get the (non-position) data into stencil order
        template<class T>
        void collectData
        (
            const GeometricField<T, fvsPatchField, surfaceMesh>& fld,
            List<List<T> >& stencilFld
        ) const;

        //- Use map to get the positions into stencil order
        template<class T>
        void collectPositions
        (
            const GeometricField<T, fvsPatchField, surfaceMesh>& fld,
            List<List<T> >& stencilFld
        ) const;
        
        //- Sum surface field contributions to create face values
        template<class Type>
        tmp<GeometricField<Vector<Type>, fvsPatchField, surfaceMesh> > weightedSum
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
            const List<List<vector> >& stencilVectorWeights
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedCentredFaceToFaceStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
