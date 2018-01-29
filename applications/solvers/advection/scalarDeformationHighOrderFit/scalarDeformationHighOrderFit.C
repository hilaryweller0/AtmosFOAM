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
    scalarDeformationHighOrderFit

Description
    Solves a transport equation for a passive scalar using explicit RK4
    time-stepping for the evolving velocity field on a plane
    for deformational flow using ghost cells to handle cyclic boundary conditions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvGhostMesh.H"
#include "deformationalFlow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    Foam::argList::validArgs.append("implicit|explicit");
    #include "setRootCase.H"
    const Switch implicit = (args.args()[1] == "implicit");
    if (!implicit && (args.args()[1] != "explicit"))
    {
        args.printUsage();
        exit(0);
    }

    #include "createTime.H"
    #include "createMesh.H"
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
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

    // Create the class for the deformational flow
    deformationalFlow defFlow
    (
        IOdictionary
        (
            IOobject
            (
                "deformationalAdvectionDict", "system", runTime,
                IOobject::MUST_READ
            )
        )
    );

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        defFlow.update(phi, U, Uf);
        defFlow.update(phiGhost);
        #include "CourantNo.H"

        for (int corr = 0; corr < nCorr; corr++)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) + 0.5*divPhiT + 0.5*divPhiT.oldTime()
            );

            if (implicit)
            {
                TEqn += 0.5*fvm::div(phi, T) - 0.5*fvc::div(phi, T);
            }

           TEqn.solve();

            // Map T to ghost mesh, calculate divergence and map back
            TGhost = ghostMesh.mapToGhost(T);
            divPhiTGhost = fvc::div(phiGhost, TGhost);
            divPhiT = ghostMesh.mapFromGhost(divPhiTGhost);
        }

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
