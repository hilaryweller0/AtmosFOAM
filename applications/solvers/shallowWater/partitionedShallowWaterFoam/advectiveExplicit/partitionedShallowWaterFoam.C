/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    partitionedShallowWaterFoamAdvExp

Description
    Transient Solver for shallow water partitioned flow - fully explicit with
    momentum equations in advective form

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #define dt runTime.deltaT()
    #include "createVariables.H"
    #include "createFields.H"
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #include "energyInit.H"
    #include "energyTransfersInit.H"
    #include "writeDiagnosticsInit.H"
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int icorr=0; icorr < nCorr; icorr++)
        {
            #include "continuityEquation.H"
            #include "updateSigma.H"
            #include "momentumEquation.H"
        }
        #include "massTransfer.H"
        #include "energy.H"
        #include "writeDiagnostics.H"
        
        runTime.write();

        Info << "sigma[0] goes from " << min(sigma[0]).value() << " to "
             << max(sigma[0]).value() << endl;
        Info << "sigma[1] goes from " << min(sigma[1]).value() << " to "
             << max(sigma[1]).value() << endl;
        Info << "h goes from " << min(h).value() << " to "
             << max(h).value() << endl;
        Info << "Total h: " << sum(h).value() << endl;
        Info << "Energy change: " << normalEnergyChange << endl;
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
