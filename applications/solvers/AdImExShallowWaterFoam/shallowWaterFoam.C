/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    AdImExShallowWaterFoam

Description
    Transient solver for inviscid shallow-water equations with rotation with
    adaptive implicit-explicit advection.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.
    
    Adaptive implicit-explicit not implemented yet. 
    
    There is an option "opSplit" to use or not use operator splitting

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvMesh.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcLaplacian.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "numericalParameters.H"
    #define dt runTime.deltaT()
    #define alpha num.alpha
    #include "readEarthProperties.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << nl << endl;

        #include "CourantNo.H"

        // Outer Corrections
        for(int outerCorr = 0; outerCorr < num.nOuterCorrs; outerCorr++)
        {
            hf = fvc::interpolate(h);
        
            // Create and solve the momentum equation
            // Rate of change of momentum with/without pressure gradient
            dhUdt = -h*(F ^ U);
            if (!num.opSplit) dhUdt -= ghGradh;

            // Explicit momentum solve
            /*dhUdt -= fvc::div(phi, U)
                  - ((fvc::div(phi, U,"div(phi,U)")) & gHat)*gHat;
            fvVectorMatrix UEqn
            (
                fvm::ddt(h,U)
             == (1-alpha)*dhUdt.oldTime() + alpha*dhUdt
            );
            UEqn.solve();*/

            // Momentum equation with implicit advection, without radial component
            fvVectorMatrix UEqn
            (
                fvm::ddt(h,U)
              + fvm::div(alpha*phi, U, "div(phi,U)")
             == (1-alpha)*dhUdt.oldTime() + alpha*dhUdt
              + ((fvc::div(alpha*phi, U,"div(phi,U)")) & gHat)*gHat
            );
            UEqn.solve();

            // Update rate of change WITHOUT pressure gradient (to be added
            // after the pressure equation)
            U -= (U & gHat)*gHat;
            dhUdt -= fvc::div(phi, U)
                  - ((fvc::div(phi, U,"div(phi,U)")) & gHat)*gHat;
            
            // Remove the pressure gradient, if it was previously included
            if (!num.opSplit) dhUdt += ghGradh;

            // The momentum without the pressure gradient
            volVectorField hU = h.oldTime() * U.oldTime()
                              + dt*((1-alpha)*dhUdt.oldTime() + alpha*dhUdt);
            // The flux without the pressure gradient
            phi = fvc::flux(hU) // Convergent version stops here
                - alpha*dt*magg*hf*fvc::snGrad(h0)*mesh.magSf();
            
            // Construct and solve the pressure equation
            for(int icorr = 0; icorr < num.nPressureCorrs; icorr++)
            {
                // Solve pressure equation
                for(int orthCorr = 0; orthCorr < num.nNonOrthogCorrs;orthCorr++)
                {
                    fvScalarMatrix hEqn
                    (
                        fvm::ddt(h)
                      + fvc::div((1-alpha)*phi.oldTime())
                      + fvc::div(alpha*phi)
                      - fvm::laplacian(sqr(alpha)*dt*magg*hf, h)
                    );
                    hEqn.solve();
                    
                    bool lastIter = icorr == num.nPressureCorrs-1 
                            && orthCorr == num.nNonOrthogCorrs-1 && alpha > 0;
                    if (lastIter) phi += hEqn.flux()/alpha;
                }
            }

            // Back substitutions
            ghGradh = fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf());
            ghGradh -= (ghGradh & gHat)*gHat;
            U = (hU - alpha*dt*ghGradh)/h;
        }
        
        dhUdt -= ghGradh;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
