/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    stratifiedBoundaryFoam

Description
    Steady-state solver for incompressible, 1D turbulent flow, typically to
    generate boundary layer conditions at an inlet, for use in a simulation.

    Boundary layer code to calculate the U, k and epsilon distributions.
    Used to create inlet boundary conditions for experimental comparisons
    for which U and k have not been measured.
    Turbulence model is runtime selectable.
    Reads in dbdz profile to set stabality functions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"
#include "wallFvPatch.H"
#include "wallDist.H"
#include "makeGraph.H"
//#include "quaternion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCaseLists.H"

    #include "createTime.H"
    const dimensionedScalar& dt = runTime.deltaT();
    #include "createMesh.H"
    #include "createFields.H"
    #include "interrogateWallPatches.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Write initial conditions
    #include "makeGraphs.H"
    makeGraph(z, dbdz, "dbdz", gFormat);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix divR(turbulence->divDevReff(U));
        
        fvVectorMatrix UEqn
        (
            fvm::ddt(U) + divR
         == (Coriolisf ^ (geostrophicWind - U)) + fvOptions(U)
        );
        UEqn.relax();
        fvOptions.constrain(UEqn);
        UEqn.solve();
        fvOptions.correct(U);

        laminarTransport.correct();
        turbulence->correct();

        #include "evaluateNearWall.H"

        if (runTime.writeTime())
        {
            #include "makeGraphs.H"
        }

        Info<< nl <<"ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
