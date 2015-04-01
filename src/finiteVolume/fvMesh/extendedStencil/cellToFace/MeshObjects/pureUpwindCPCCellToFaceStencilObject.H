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
    Foam::pureUpwindCPCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pureUpwindCPCCellToFaceStencilObject_H
#define pureUpwindCPCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CPCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pureUpwindCPCCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class pureUpwindCPCCellToFaceStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        pureUpwindCPCCellToFaceStencilObject
    >,
    public extendedUpwindCellToFaceStencil
{

public:

    TypeName("pureUpwindCPCCellToFaceStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit pureUpwindCPCCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                pureUpwindCPCCellToFaceStencilObject
            >(mesh),
            extendedUpwindCellToFaceStencil(CPCCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated pure upwind stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


    //- Destructor
    virtual ~pureUpwindCPCCellToFaceStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
