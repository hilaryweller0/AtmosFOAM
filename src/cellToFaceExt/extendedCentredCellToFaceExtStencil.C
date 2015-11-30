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

\*---------------------------------------------------------------------------*/

#include "mapDistribute.H"
#include "extendedCentredCellToFaceExtStencil.H"
#include "cellToFaceExtStencil.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCentredCellToFaceExtStencil::extendedCentredCellToFaceExtStencil
(
    const cellToFaceExtStencil& stencil
)
:
    extendedCellToFaceExtStencil(stencil.mesh()),
    untransformedElements_(stencil.untransformedElements())
{
    // Calculate distribute map (also renumbers elements in stencil)
    List<Map<label> > compactMap(Pstream::nProcs());
    mapPtr_.reset
    (
        new mapDistribute
        (
            stencil.globalNumbering(),
            untransformedElements_,
            stencil.mesh().globalData().globalTransforms(),
            stencil.transformedElements(),
            transformedElements_,
            compactMap
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Per face which elements of the stencil to keep.
void Foam::extendedCentredCellToFaceExtStencil::compact()
{
    boolList isInStencil(map().constructSize(), false);

    forAll(untransformedElements_, faceI)
    {
        // Mark untransformed stencil elements
        {
            const labelList& stencilCells = untransformedElements_[faceI];

            forAll(stencilCells, i)
            {
                isInStencil[stencilCells[i]] = true;
            }
        }

        // Mark transformed stencil elements
        {
            const labelList& transformedElements = transformedElements_[faceI];
    
            forAll(transformedElements, i)
            {
                isInStencil[transformedElements[i]] = true;
            }
        }
    }

    mapPtr_().compact(isInStencil, Pstream::msgType());
}


// ************************************************************************* //
