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
    Foam::upwindCPCCellToFaceStencilNewObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCPCCellToFaceStencilNewObject_H
#define upwindCPCCellToFaceStencilNewObject_H

#include "extendedUpwindCellToFaceStencilNew.H"
#include "CPCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class upwindCPCCellToFaceStencilNewObject Declaration
\*---------------------------------------------------------------------------*/

class upwindCPCCellToFaceStencilNewObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        upwindCPCCellToFaceStencilNewObject
    >,
    public extendedUpwindCellToFaceStencilNew
{

public:

    TypeName("upwindCPCCellToFaceStencilNew");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit upwindCPCCellToFaceStencilNewObject
        (
            const fvMesh& mesh,
            const bool pureUpwind,
            const scalar minOpposedness,
            const bool linearCorr = true
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                upwindCPCCellToFaceStencilNewObject
            >(mesh),
            extendedUpwindCellToFaceStencilNew
            (
                CPCCellToFaceStencil(mesh),
                pureUpwind,
                minOpposedness,
                linearCorr
            )
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated off-centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


    //- Destructor
    virtual ~upwindCPCCellToFaceStencilNewObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
