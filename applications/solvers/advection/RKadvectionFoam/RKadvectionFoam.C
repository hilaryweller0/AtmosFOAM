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
    RKadvectionFoam

Description
    Solve the linear advection equation for constant in time wind field with
    an explicit RK scheme

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "butcherTableau.H"
#include "RKfield.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "createRKFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
        
        for (int iRK = 0; iRK < Bt.nSteps(); iRK++)
        {
            // Explicit RK stage
            T = T.oldTime() + dt*Bt.RKstep(iRK, dTdt);

            // update time derivatives for next stage
            dTdt[iRK] = -fvc::div(phi, T);
        }

        if (max(mag(Bt.weights())) > SMALL)
        {
            // Final RK steps
            T = T.oldTime() + dt*Bt.RKfinal(dTdt);
        }

        Info << "T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << nl << endl;
        runTime.write();
                
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
