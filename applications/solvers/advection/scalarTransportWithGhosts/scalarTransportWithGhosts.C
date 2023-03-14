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
    scalarTransportWithGhosts

Description
    Solves a transport equation for a passive scalar using Implicit
    time-stepping on a periodic domain with ghost cells

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvGhostMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solution().subDict("iterations");
    const int nCorr = readLabel(itsDict.lookup("nCorr"));
    
    // Create the ghost mesh
    fvGhostMesh ghostMesh
    (
        IOobject
        (
            "ghostMesh", runTime.timeName(), runTime, IOobject::MUST_READ
        ),
        mesh
    );
    
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;
    
    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        #include "CourantNo.H"

        for (int corr = 0; corr < nCorr; corr++)
        {
            solve
            (
                fvm::ddt(T)
              + 0.5*fvm::div(phi, T)
              - 0.5*fvc::div(phi, T)
              + 0.5*divPhiT
              + 0.5*divPhiT.oldTime()
            );

            // Map T to ghost mesh, calculate divergence and map back
            TGhost = ghostMesh.mapToGhost(T);
            divPhiTGhost = fvc::div(phiGhost, TGhost);
            divPhiT = ghostMesh.mapFromGhost(divPhiTGhost);
        }

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
