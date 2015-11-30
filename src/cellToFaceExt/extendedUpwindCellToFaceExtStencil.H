/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2015-11-17 AtmosFOAM, Hilary Weller, University of Reading added support for
 cyclic boundaries
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
    Foam::extendedUpwindCellToFaceExtStencil

Description
    Creates upwind stencil by shifting a centred stencil to upwind and downwind
    faces and optionally removing all non-(up/down)wind faces ('pureUpwind').

    Note: the minOpposedness parameter is to decide which upwind and
    downwind faces to combine the stencils from. If myArea is the
    local area and upwindArea
    the area of the possible upwind candidate it will be included if
        (upwindArea & myArea)/magSqr(myArea) > minOpposedness
    so this includes both cosine and area. WIP.

SourceFiles
    extendedUpwindCellToFaceExtStencil.C
    extendedUpwindCellToFaceExtStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedUpwindCellToFaceExtStencil_H
#define extendedUpwindCellToFaceExtStencil_H

#include "extendedCellToFaceExtStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToFaceExtStencil;
class mapDistribute;

/*---------------------------------------------------------------------------*\
               Class extendedUpwindCellToFaceExtStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedUpwindCellToFaceExtStencil
:
    public extendedCellToFaceExtStencil
{
    // Private data

        //- Does stencil contain upwind points only
        const bool pureUpwind_;

        //- Swap map for getting neigbouring data
        autoPtr<mapDistribute> ownMapPtr_;
        autoPtr<mapDistribute> neiMapPtr_;

        //- Per face the owner and neighbour stencils.
        labelListList ownStencil_;
        labelListList neiStencil_;

        //- Per face the owner and neighbour transformed elements
        List<labelPairList> ownTransformedElements_;
        List<labelPairList> neiTransformedElements_;

    // Private Member Functions

        //- Find most 'opposite' faces of cell
        void selectOppositeFaces
        (
            const boolList& nonEmptyFace,
            const scalar minOpposedness,
            const label faceI,
            const label cellI,
            DynamicList<label>& oppositeFaces
        ) const;

        //- Transport (centred) face stencil to 'opposite' face.
        void transportStencil
        (
            const boolList& nonEmptyFace,
            const labelListList& faceStencil,
            const List<labelPairList>& transformedFaceStencil,
            const scalar minOpposedness,
            const label faceI,
            const label cellI,
            const bool stencilHasNeighbour,

            DynamicList<label>& oppositeFaces,
            labelHashSet& faceStencilSet,
            HashSet<labelPair>& transformedFaceStencilSet,
            labelList& transportedStencil,
            labelPairList& transportedTransformedStencil
        ) const;

        //- Transport (centred) face stencil to 'opposite' faces.
        void transportStencils
        (
            const labelListList& faceStencil,
            const List<labelPairList>& transformedStencil,
            const scalar minOpposedness,
            labelListList& ownStencil,
            labelListList& neiStencil,
            List<labelPairList>& ownTransformedElements,
            List<labelPairList>& neiTransformedElements
        );


        //- Disallow default bitwise copy construct
        extendedUpwindCellToFaceExtStencil(const extendedUpwindCellToFaceExtStencil&);

        //- Disallow default bitwise assignment
        void operator=(const extendedUpwindCellToFaceExtStencil&);


public:

    // Constructors

        //- Construct from mesh and uncompacted centred face stencil.
        //  Transports facestencil to create owner and neighbour versions.
        //  pureUpwind to remove any remaining downwind cells.
        extendedUpwindCellToFaceExtStencil
        (
            const cellToFaceExtStencil&,
            const bool pureUpwind,
            const scalar minOpposedness
        );


    // Member Functions

        bool pureUpwind() const
        {
            return pureUpwind_;
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& ownMap() const
        {
            return ownMapPtr_();
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& neiMap() const
        {
            return neiMapPtr_();
        }

        //- Return reference to the stencil
        const labelListList& ownStencil() const
        {
            return ownStencil_;
        }

        //- Return reference to the stencil
        const labelListList& neiStencil() const
        {
            return neiStencil_;
        }

        //- Return reference to the transformed owner stencil
        const List<labelPairList>& ownTransformedElements() const
        {
            return ownTransformedElements_;
        }

        //- Return reference to the transformed neighbour stencil
        const List<labelPairList>& neiTransformedElements() const
        {
            return neiTransformedElements_;
        }

        //- Use map to get the positions into stencil order for owner
        template<class T>
        void collectOwnPositions
        (
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T> >& stencilFld
        ) const
        {
            extendedCellToFaceExtStencil::collectData
            (
                ownMap(),
                ownStencil(),
                ownTransformedElements(),
                fld(),
                stencilFld(),
                mapDistribute::transformPosition()
            );
        }

        //- Use map to get the positions into stencil order for neighbour
        template<class T>
        void collectNeiPositions
        (
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T> >& stencilFld
        ) const
        {
            extendedCellToFaceExtStencil::collectData
            (
                neiMap(),
                neiStencil(),
                neiTransformedElements(),
                fld(),
                stencilFld(),
                mapDistribute::transformPosition()
            );
        }

        //- Sum vol field contributions to create face values
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > weightedSum
        (
            const surfaceScalarField& phi,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar> >& ownWeights,
            const List<List<scalar> >& neiWeights
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedUpwindCellToFaceExtStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
