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

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int icorr=0; icorr < nCorr; icorr++)
        {
            // Advect h in each partition
            for(label ip = 0; ip < nParts; ip++)
            {
                h[ip] = h[ip].oldTime() - dt*fvc::div(volFlux[ip], h[ip]);
            
                if (ip == 0) hSum = h[ip];
                else hSum += h[ip];
            }
            // Update sigma (diagnostic)
            for(label ip = 0; ip < nParts; ip++)
            {
                sigma[ip] = h[ip]/hSum;
            }
            
            // Update the velocity in each partition
            surfaceScalarField ggradh = g*fvc::snGrad(hSum)*mesh.magSf();
            for (int ucorr = 0; ucorr < nUCorr; ucorr++)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    volFlux[ip] = volFlux[ip].oldTime() - dt*
                    (
                        ((Uf[ip]&fvc::interpolate(fvc::grad(Uf[ip])))&mesh.Sf())
                      + ((twoOmegaf^Uf[ip]) & mesh.Sf())
                      + ggradh
                    );
                
                    u[ip] = fvc::reconstruct(volFlux[ip]);
                    Uf[ip] = fvc::interpolate(u[ip]);
                    Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                            *mesh.Sf()/sqr(mesh.magSf());
                }
            }
        }

        #include "energy.H"
        Info << "sigma[0] goes from " << min(sigma[0]).value() << " to "
             << max(sigma[0]).value() << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
