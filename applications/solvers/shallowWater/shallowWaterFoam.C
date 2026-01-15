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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvMesh.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
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
            // Create and solve the momentum equation
            // Rate of change without implicit advection
            dhUdt = -h*(F ^ U) - magg*h*fvc::grad(h + h0);
            
            // How close is the initial gradient wind balance
            volScalarField gradWind("gradWind", mag(fvc::div(phi,U) + h*(F^U)));
            gradWind.write();
            volScalarField ghGradh("ghGradh", mag(magg*h*fvc::grad(h + h0)));
            ghGradh.write();
            volScalarField imbalance("imbalance", gradWind-ghGradh);
            imbalance.write();
            
            // Momentum equation with implicit advection
            fvVectorMatrix UEqn
            (
                fvm::ddt(h,U)
              + fvm::div(alpha*phi, U, "div(phi,U)")
             == (1-alpha)*dhUdt.oldTime() + alpha*dhUdt
            );
            UEqn.solve();

            // Update rate of change without pressure gradient (to be added
            // after the pressure equation)
            dhUdt -= fvc::div(phi, U) + magg*h*fvc::grad(h+h0);

            // Constrain the momentum to be in the geometry if 3D geometry
            if (mesh.nGeometricD() == 3)
            {
                U -= (gHat & U)*gHat;
                dhUdt -= (gHat & dhUdt)*gHat;
                U.correctBoundaryConditions();
                dhUdt.correctBoundaryConditions();
            }
            
            // Calculate the momentum without the pressure gradient
            volVectorField hU = h.oldTime() * U.oldTime()
                              + dt*((1-alpha)*dhUdt.oldTime() + alpha*dhUdt);
            
            // Convert the momentum into a flux
            phi = fvc::flux(hU);
            
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
                      - fvm::laplacian(alpha*dt*magg*h, h)
                      - fvc::laplacian(alpha*dt*magg*h, h0)
                    );
                    hEqn.solve();

                    // Back substitutions
                    if (orthCorr == num.nNonOrthogCorrs-1
                        && icorr == num.nPressureCorrs-1)
                    {
                        surfaceScalarField hFlux = hEqn.flux()
                            - alpha*dt*magg*fvc::interpolate(h)*fvc::snGrad(h0)*mesh.magSf();
                        phi += hFlux;
                        volVectorField hUinc = fvc::reconstruct(hFlux);
                        hU  += hUinc;
                        dhUdt += hUinc/(alpha*dt);
                        dhUdt.correctBoundaryConditions();
                    }
                }

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                }

                U = hU/h;
                U.correctBoundaryConditions();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
