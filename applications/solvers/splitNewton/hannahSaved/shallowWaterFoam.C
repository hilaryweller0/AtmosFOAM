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
    ColocShallowWaterFoamSplitNewton

Description
    Transient solver for inviscid shallow-water equations with rotation with
    an RK IMEX scheme.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.
    

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
    #include "readEarthProperties.H"
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalar totalE = gSum((0.5*h*magSqr(U) + 0.5*magg*sqr(h))().primitiveField() * mesh.V());
    Info<< "E = " << totalE << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << nl << endl;

        #include "CourantNo.H"
        
        // half gravity (expl)
                
        h = h.oldTime() - 0.5*dt*fvc::div(phi); 
        hf = fvc::interpolate(h);
        
        volVectorField ghGradh = fvc::reconstruct(magg*hf*fvc::snGrad(h.oldTime()+h0)*mesh.magSf());
        ghGradh -= (ghGradh & gHat) * gHat; // outer velocity correction
        
        U = 
        (
            h.oldTime()*U
          - 0.5*dt*ghGradh
        )/h;
        
        
        // half Coriolis
        
        U.oldTimeRef() = U;
        
        for(int it = 0; it < 2; it++)
        {
            U = U.oldTime() - 0.5*dt*(F ^ U);
        }
        
        phi = fvc::flux(h*U);
        
        
        // Momentum advection with iterations for the non-linearity
        
        // First store the previous velocity (from after half a Coriolis + gravity)
        U.oldTimeRef() = U;
        
        // First implicit RK stage
        for(int its = 0; its < num.nItsU1; its++)
        {
            fvVectorMatrix UEqn
            (
                1/dt* (fvm::Sp(h,U_prime) + h*U - h*U.oldTime())  //1/dt * (fvm::Sp(h,U) - h*U.oldTime())
              + 0.25*fvc::div(phi, U)
              + 0.5*fvm::div(phi, U_prime)
              - 0.25*(fvc::div(phi, U) & gHat)*gHat // outer velocity correction
              - 0.5*(fvc::div(phi, U_prime) & gHat)*gHat // outer velocity correction
            );
            UEqn.solve();
            U += U_prime;
            phi = fvc::flux(h*U);
        };
        
        volVectorField dhUdt_a = fvc::div(phi, U, "div(phi,U)");
        dhUdt_a -= (dhUdt_a & gHat)*gHat;
        
        // Second implicit RK stage
        for(int its = 0; its < num.nItsU2; its++)
        {
            fvVectorMatrix UEqn
            (
                1/dt* (fvm::Sp(h,U_prime) + h*U - h*U.oldTime()) //1/dt * (fvm::Sp(h,U) - h*U.oldTime())
              + 2./3.* dhUdt_a
              + 1./3.*fvc::div(phi, U)
              + 2./3.*fvm::div(phi, U_prime)
              - 1./3.*(fvc::div(phi, U)& gHat)*gHat // outer velocity correction
              - 2./3.*(fvc::div(phi, U_prime) & gHat)*gHat // outer velocity correction
            );
            UEqn.solve();
            U += U_prime;
            phi = fvc::flux(h*U);
        };
        
        
        // Final half Coriolis
        
        U -= 0.5*dt*(F ^ U);
        
        
        // Previous stages are no longer needed
        U.oldTimeRef() = U;
        h.oldTimeRef() = h;
        phi = fvc::flux(h*U);
        
        
        // Final gravity (implicit), with non-linear iterations
        
        for(int its = 0; its < num.nItsh; its++)
        {
            hf = fvc::interpolate(h);
            fvScalarMatrix hEqn
            (
                1/dt* (fvm::Sp(1,h_prime) + h - h.oldTime()) //1/dt * (fvm::Sp(1,h) - h.oldTime())
              + 0.5*fvc::div(phi)
              + 0.5*fvm::div(linearInterpolate(U) & mesh.Sf(),h_prime,"div(phi,h)")
              - 0.25*fvc::laplacian(dt*magg*hf, h)
              - 0.25*fvm::laplacian(dt*magg*hf, h_prime)
              - 0.25*fvm::div(dt*magg*fvc::snGrad(h)*mesh.magSf(), h_prime,"div(phi,h)")
              - 0.25*fvc::laplacian(dt*magg*hf, h0)
            );
            hEqn.solve();
            h += h_prime;
            //h_prime = 0;
        };
        
        ghGradh = fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf());
        ghGradh -= (ghGradh & gHat) * gHat; // outer velocity correction
        
        U = 
        (
            h.oldTime()*U.oldTime() 
          - 0.5 * dt * ghGradh
        )/h;
        phi = fvc::flux(h*U);
        
        scalar totalE = gSum((0.5*h*magSqr(U) + 0.5*magg*sqr(h))().primitiveField() * mesh.V());
        Info<< "E = " << totalE << endl;
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
