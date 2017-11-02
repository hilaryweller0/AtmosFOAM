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
    testTetDecompose

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    forAll(mesh.C(), cellI)
    {
        List<tetIndices> tets = polyMeshTetDecomposition::cellTetIndices(mesh, cellI);
        forAll(tets, tetI)
        {
            tetPointRef tet = tets[tetI].tet(mesh);
            for (label tetFaceI = 0; tetFaceI < 4; tetFaceI++)
            {
                triPointRef tri = tet.tri(tetFaceI);
                Info << tetI << " " << tri.a().x() << " " << tri.a().y() << " " << tri.a().z() << endl;
                Info << tetI << " " << tri.b().x() << " " << tri.b().y() << " " << tri.b().z() << endl;
                Info << tetI << " " << tri.c().x() << " " << tri.c().y() << " " << tri.c().z() << endl;
                Info << tetI << " " << tri.a().x() << " " << tri.a().y() << " " << tri.a().z() << endl;
                Info << endl;
            }
            Info << endl;
        }
        break;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
