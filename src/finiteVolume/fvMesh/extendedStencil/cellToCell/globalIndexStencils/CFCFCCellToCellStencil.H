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

Class
    Foam::CFCFCCellToCellStencil

Description
    Cell stencil around each cell containing all face neighbours of face
    neighbours

SourceFiles
    CFCFCCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef CFCFCCellToCellStencil_H
#define CFCFCCellToCellStencil_H

#include "cellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CFCFCCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class CFCFCCellToCellStencil
:
    public cellToCellStencil
{
    // Private Member Functions

        void calcFaceBoundaryData(labelList& neiGlobal) const;

        void calcCellStencil(labelListList& globalCellCells) const;


public:

    // Constructors

        //- Construct from mesh
        explicit CFCFCCellToCellStencil(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
