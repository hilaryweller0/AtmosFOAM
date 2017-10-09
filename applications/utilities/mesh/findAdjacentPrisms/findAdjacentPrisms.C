/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

Application
    findAdjacentPrisms

Description
    Write a faceSet containing interior faces that connect two prisms.
    Designed to be used in conjunction with slantedCellMesh to improve
    mesh quality by merging triangular (prismatic) slanted cells found
    near the lower boundary of slanted cell meshes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faceSet.H"
#include "prismMatcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    faceSet adjacentPrismSet(mesh, "adjacentPrismSet", IOobject::NO_READ, IOobject::AUTO_WRITE);
    prismMatcher prism;

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        const label& own = mesh.faceOwner()[faceI];
        const label& neighbour = mesh.faceNeighbour()[faceI];

        if (prism.isA(mesh, own) && prism.isA(mesh, neighbour))
        {
            adjacentPrismSet.insert(faceI);
        }
    }

    adjacentPrismSet.write();

    Info << "Written " << adjacentPrismSet.size() << " faces with two adjacent prisms to adjacentPrismSet" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
