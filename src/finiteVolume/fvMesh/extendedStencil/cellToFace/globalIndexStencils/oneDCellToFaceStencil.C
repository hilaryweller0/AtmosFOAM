/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "oneDCellToFaceStencil.H"
#include "CFCCellToCellStencil.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDCellToFaceStencil::oneDCellToFaceStencil(const polyMesh& mesh)
:
    cellToFaceStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    CFCCellToCellStencil globalCellCells(mesh);

    // Add stencils of neighbouring cells to create faceStencil
    labelListList faceStencil;
    calcFaceStencil(globalCellCells, faceStencil);
    
    // Remove cells that are not in line with the face
    forAll(faceStencil, faceI)
    {
        // Position and normal vector of this face
        const point& Cf = mesh.faceCentres()[faceI];
        const vector& Sf = mesh.faceAreas()[faceI];
        const label oldSize = faceStencil[faceI].size();
        label newSize = 0;
        labelList inlineStencil(oldSize);
        forAll(faceStencil[faceI], ic)
        {
            label cellI = faceStencil[faceI][ic];
            // Only consider cell, not boundary faces
            if (cellI < mesh.nCells())
            {
                const point& C = mesh.cellCentres()[cellI];
                if ( mag(mag((C-Cf) & Sf)/(mag(C-Cf)*mag(Sf)) - 1) < SMALL)
                {
                    inlineStencil[newSize++] = cellI;
                }
            }
        }
        
        if (newSize == 0 && faceI < mesh.nInternalFaces())
        {
            FatalErrorIn("oneDCellToFaceStencil(mesh)") << " face " << faceI
                << " has no stencil in line" << exit(FatalError);
        }
        
        inlineStencil.setSize(newSize);
        faceStencil[faceI] = inlineStencil;
    }

    // Transfer to *this
    transfer(faceStencil);
}


// ************************************************************************* //
