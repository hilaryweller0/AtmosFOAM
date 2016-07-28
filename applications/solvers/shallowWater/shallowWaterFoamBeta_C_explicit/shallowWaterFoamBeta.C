/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    shallowWaterFoamBeta_C_explicit

Description
    Transient explicit solver for inviscid shallow-water equations in
    flux form with rotation. Uses the beta plane approximation and a C-grid.
    The geometry needs to be in the x-y plane as rotation is assumed about
    the z-axis. Gravity waves and advection treated explicitly

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
    // Read solution controls
    const int nCorr = readLabel(mesh.solutionDict().lookup("nCorrectors"));
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "CourantNo.H"

        // --- Outer loop

        for (int corr=0; corr<nCorr; corr++)
        {
            // Calculate the explicit source terms in the momentum equations
            // without Coriolis
            dFluxdt = -hf*gsnGradh0
                    - (fvc::interpolate(fvc::div(flux,U)) & mesh.Sf());
            
            // Crank Nicolson update of flux without implicit terms and
            // without Coriolis which needs to be balanced at the boundary
            fluxU = flux.oldTime() + runTime.deltaT()*
            (
                (1-offCentre)*dFluxdt.oldTime() + offCentre*dFluxdt
            );
            
            // Add components of the flux that can create flow over boundaries
            flux = fluxU - dt*offCentre*hf*g*fvc::snGrad(h)*mesh.magSf();
            if (rotating) flux -= dt*offCentre*hf*(fSfxk & Uf);

            // Solve the continuity equation
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvc::div((1-offCentre)*flux.oldTime() + offCentre*flux)
            );

            hEqn.solve();

            // update diagnsotic fields
            hf = fvc::interpolate(h);
            U = fvc::reconstruct(flux/hf);
            Uf = fvc::interpolate(U);
        }
        // Rate of change of flux for next time-step (including all terms)
        dFluxdt -= g*hf*fvc::snGrad(h)*mesh.magSf();
        if (rotating) dFluxdt -= hf*(fSfxk & Uf);

        //hTotal == h + h0;
        runTime.write();

        Info<< "\n    ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
            << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
