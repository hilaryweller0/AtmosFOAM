/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"
#include "SubList.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "OSspecific.H"
#include "demandDrivenData.H"

#include "pointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::calcDirections() const
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        solutionD_[cmpt] = 1;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::polyMesh::nGeometricD() const
{
    label nD = 3;
    label nEmptyPatches = 0;

    forAll(boundaryMesh(), patchi)
    {
        if (boundaryMesh()[patchi].size())
        {
            if (isA<emptyPolyPatch>(boundaryMesh()[patchi]))
            {
                nEmptyPatches++;
            }
        }
    }

    nD -= label(nEmptyPatches/2);

    return nD;
}


// ************************************************************************* //
