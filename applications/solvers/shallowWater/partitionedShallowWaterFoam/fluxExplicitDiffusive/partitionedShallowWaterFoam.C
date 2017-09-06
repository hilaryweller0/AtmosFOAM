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
    partitionedShallowWaterFoamFluxExp

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
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 2);
    const int nUCorr = itsDict.lookupOrDefault<int>("nUCorrs", 2);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #include "energyInit.H"
    #include "writeDiagnosticsInit.H"
    
    const dimensionedScalar K("K",dimensionSet(0,2,-1,0,0),scalar(5000000));

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int icorr=0; icorr < nCorr; icorr++)
        {
            // Advect h in each partition
            for(label ip = 0; ip < nParts; ip++)
            {
                h[ip] = h[ip].oldTime() - dt*
                (
                    fvc::div(flux[ip])
                );
                for(label ip2 = 0; ip2 < nParts; ip2++)
                {
                    h[ip] -= (2*(ip2 != ip)-1)*dt*K*hOld[ip2]*fvc::laplacian(sigma[ip2]);
                }
                
                hf[ip] = fvc::interpolate(h[ip]);
                
                Info << "h[" << ip << "] goes from " 
                     << min(h[ip].internalField()).value() 
                     << " to " 
                     << max(h[ip].internalField()).value() << endl;

                if (ip == 0) hSum = h[ip];
                else hSum += h[ip];
            }
            // Update sigma (diagnostic)
            for(label ip = 0; ip < nParts; ip++)
            {
                sigma[ip] = h[ip]/hSum;
                hOld[ip] = h[ip];
            }
            
            // Update the velocity in each partition
            surfaceScalarField ggradh = g*fvc::snGrad(hSum)*mesh.magSf();
            for (int ucorr = 0; ucorr < nUCorr; ucorr++)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    flux[ip] = flux[ip].oldTime() - dt*
                    (
                        fvc::flux(fvc::div(flux[ip], u[ip]))
                      //+ hf[ip]*((twoOmegaf^Uf[ip]) & mesh.Sf())
                      + hf[ip]*ggradh
                    );
                    
                    for(label ip2 = 0; ip2 < nParts; ip2++)
                    {
                        flux[ip] -= (2*(ip2 != ip)-1)*dt*fvc::interpolate(K*hOld[ip2]*fvc::laplacian(sigma[ip2]))*(Uf[ip2] & mesh.Sf());
                    }
                    
                    u[ip] = fvc::reconstruct(flux[ip]/hf[ip]);
                    //u[ip] = fvc::reconstruct(flux[ip])/h[ip];
                    Uf[ip] = fvc::interpolate(u[ip]);
                    Uf[ip] += (flux[ip]/hf[ip] - (Uf[ip] & mesh.Sf()))
                            *mesh.Sf()/sqr(mesh.magSf());
                }
            }
        }

        #include "energy.H"
        #include "writeDiagnostics.H"
        
        Info << "sigma[0] goes from " << min(sigma[0]).value() << " to "
             << max(sigma[0]).value() << endl;
        Info << "sigma[1] goes from " << min(sigma[1]).value() << " to "
             << max(sigma[1]).value() << endl;
        Info << "h goes from " 
             << min(h[0].internalField()+h[1].internalField()).value() 
             << " to " 
             << max(h[0].internalField()+h[1].internalField()).value() 
             << endl;
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
