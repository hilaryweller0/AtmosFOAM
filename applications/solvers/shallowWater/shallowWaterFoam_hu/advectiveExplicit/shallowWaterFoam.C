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
    shallowWaterFoamAdvExp

Description
    Transient Solver for shallow water flow - fully explicit with
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
    
    const dictionary& itsDict = mesh.solution().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nUCorr = itsDict.lookupOrDefault<int>("nUCorrs", 1);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #include "energyInit.H"
    #include "writeDiagnosticsInit.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "courantNo.H"

        for (int icorr=0; icorr < nCorr; icorr++)
        {
            h = h.oldTime() - dt*
            (
                fvc::div(volFlux, h)
            );
            
            // Update the velocity in each partition
            for (int ucorr = 0; ucorr < nUCorr; ucorr++)
            {
                surfaceScalarField ggradh = g*fvc::snGrad(h)*mesh.magSf();
                
                //Update prognostic variables.
                volFlux = volFlux.oldTime() - dt*
                (
                    ((Uf&fvc::interpolate(fvc::grad(Uf)))&mesh.Sf())
                    
                  //+ ((twoOmegaf^Uf) & mesh.Sf())
                  + ggradh
                );
                
                U = fvc::reconstruct(volFlux);
                Uf = fvc::interpolate(U);
                Uf += (volFlux - (Uf & mesh.Sf()))
                        *mesh.Sf()/sqr(mesh.magSf());
            }
        }

        #include "energy.H"
        #include "writeDiagnostics.H"
        Info << "h goes from " << min(h).value() << " to "
             << max(h).value() << endl;
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
