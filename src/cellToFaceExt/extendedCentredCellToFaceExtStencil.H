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
    Foam::extendedCentredCellToFaceExtStencil

Description

SourceFiles
    extendedCentredCellToFaceExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToFaceExtStencil_H
#define extendedCentredCellToFaceExtStencil_H

#include "extendedCellToFaceExtStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToFaceExtStencil;
class mapDistribute;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToFaceExtStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToFaceExtStencil
:
    public extendedCellToFaceExtStencil
{
    // Private data

        //- Per face the stencil.
        labelListList untransformedElements_;

        //- Per face the transformed elements
        labelListList transformedElements_;

        //- Swap map for getting neigbouring data
        autoPtr<mapDistribute> mapPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedCentredCellToFaceExtStencil
        (
            const extendedCentredCellToFaceExtStencil&
        );

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToFaceExtStencil&);


public:

    // Constructors

        //- Construct from uncompacted face stencil
        explicit extendedCentredCellToFaceExtStencil(const cellToFaceExtStencil&);


    // Member Functions

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            return mapPtr_();
        }

        //- Return reference to the stencil
        const labelListList& untransformedElements() const
        {
            return untransformedElements_;
        }
        
        //- Return reference to the transformed elements in the stencil
        const labelListList& transformedElements() const
        {
            return transformedElements_;
        }

        //- After removing elements from the stencil adapt the schedule (map).
        void compact();

        //- Use map to get the data into stencil order
        template<class T>
        void collectData
        (
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T> >& stencilFld
        ) const
        {
            extendedCellToFaceExtStencil::collectData
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld,
                stencilFld,
                mapDistribute::transform()
            );
        }
        
        //- Use map to get the positions into stencil order
        template<class T>
        void collectPositions
        (
            const GeometricField<T, fvPatchField, volMesh>& fld,
            List<List<T> >& stencilFld
        ) const
        {
            extendedCellToFaceExtStencil::collectData
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld(),
                stencilFld(),
                mapDistribute::transformPosition()
            );
        }

        //- Sum vol field contributions to create face values
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar> >& stencilWeights
        ) const
        {
            return extendedCellToFaceExtStencil::weightedSum
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld,
                stencilWeights,
                mapDistribute::transform()
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
