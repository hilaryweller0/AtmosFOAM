/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "CPCCellToFaceExtStencil.H"
#include "CPCCellToCellExtStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToFaceExtStencil::CPCCellToFaceExtStencil(const polyMesh& mesh)
:
    cellToFaceExtStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    CPCCellToCellExtStencil globalCellCells(mesh);

    // Add stencils of neighbouring cells to create faceStencil
    calcFaceStencil
    (
        globalCellCells.untransformedElements(),
        globalCellCells.transformedElements(),

        untransformedElements_,
        transformedElements_
    );
}


// ************************************************************************* //
