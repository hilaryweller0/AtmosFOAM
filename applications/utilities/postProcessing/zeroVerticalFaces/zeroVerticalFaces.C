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

Application
    zeroVerticalFaces

Description
    Sets values on vertical Tf faces to zero.  Useful for comparing results
    on Charney--Phillips advection tests.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gravity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    Foam::argList::validArgs.append("fieldName");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    gravity g(mesh);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        surfaceScalarField Tf
        (
            IOobject
            (
                args.args()[1],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info  << "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            bool horizontal = mag(g.unitFaceNormal()[faceI]) > 1e-12;
            if (!horizontal) Tf[faceI] = 0;
        }
        Tf.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
