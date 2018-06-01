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
    partitionedShallowWaterFoamFluxAgrid

Description
    Transient Solver for shallow water partitioned flow - fully explicit with
    solution for velocity and h in each partition using flux form.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "createVariables.H"
    
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
            // Advect h in each partition
            for(label ip = 0; ip < nParts; ip++)
            {
                sigmah[ip] = sigmah[ip].oldTime() - dt*
                (
                    (1-offCentre)*fvc::div(volFlux[ip].oldTime(),sigmah[ip].oldTime())
                  + offCentre*fvc::div(volFlux[ip],sigmah[ip])
                );
                
                sigmahf[ip] = fvc::interpolate(sigmah[ip],"linear");
                
                Info << "sigmah[" << ip << "] goes from " 
                     << min(sigmah[ip].internalField()).value() 
                     << " to " 
                     << max(sigmah[ip].internalField()).value() << endl;
                
                if (ip == 0) h = sigmah[ip];
                else h += sigmah[ip];
            
                hu[ip] = sigmah[ip]*u[ip];
            }
            
            // Update sigma (diagnostic)
            for(label ip = 0; ip < nParts; ip++)
            {
                sigma[ip] = sigmah[ip]/h;
            }
            
            // Update the velocity in each partition
            
            //volVectorField ggradh = g*fvc::reconstruct(fvc::snGrad(hSum)*mesh.magSf());
            //volVectorField ggradhOld = g*fvc::reconstruct(fvc::snGrad(hSum)*mesh.magSf());
            for (int ucorr = 0; ucorr < nUCorr; ucorr++)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    volVectorField ggradh = g*fvc::grad(h,partName[ip]);
                    volVectorField ggradhOld = g*fvc::grad(h.oldTime(),partName[ip]);
                    hu[ip] = hu[ip].oldTime() - dt*
                    (
                        (1-offCentre)*fvc::div(volFlux[ip].oldTime(),hu[ip].oldTime())
                      + offCentre*fvc::div(volFlux[ip],hu[ip])
                      + (1-offCentre)*sigmah[ip].oldTime()*ggradhOld
                      + offCentre*sigmah[ip]*ggradh
                    );
                    
                    if (useBuoyancy) 
                    { 
                        hu[ip] += dt*sigmah[ip]*buoyancyMagnitude*momentumSource[ip]*yNorm;
                    }
                    
                    u[ip] = hu[ip]/sigmah[ip];
                    
                    Uf[ip] = fvc::interpolate(u[ip]);
                    
                    volFlux[ip] = Uf[ip] & mesh.Sf();
                    
                    flux[ip] = fvc::interpolate(hu[ip]) & mesh.Sf();
                }
            }
            
            #include "massTransfer.H"
        }

        #include "energy.H"
        #include "writeDiagnostics.H"
        
        
        Info << "Total h: " << sum(h).value() << endl;
        Info << "Energy change: " 
             << normalEnergyChange << endl;
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
