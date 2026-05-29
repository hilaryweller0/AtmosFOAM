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
    CgridShallowWaterFoam

Description
    Semi-implicit solver for inviscid non-linear shallow-water equations
    with rotation on a C-grid.

    If the geometry is 3D then it is assumed to be one layers of cells and
    the component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"

#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcReconstruct.H"
#include "fvcLaplacian.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "1",
        "solve using the first splitting scheme"
    );
    argList::addBoolOption
    (
        "2",
        "solve using the second splitting scheme"
    );
    argList::addBoolOption
    (
        "3",
        "solve using the third splitting scheme"
    );
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "readEarthProperties.H"
    #include "createFields.H"

    const int nIters = readLabel(mesh.solution().lookup("nIterations"));
    const scalar alpha =readScalar(mesh.solution().lookup("timeOffCentre"));
    
    const bool scheme1 = args.optionFound("1");
    const bool scheme2 = args.optionFound("2");
    const bool scheme3 = args.optionFound("3");
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << endl;

        #include "CourantNo.H"
        
        volScalarField h_n = h.oldTime();
        volVectorField U_n = U.oldTime();

        // Outer Iterations
        for (int iIt=0; iIt < nIters; iIt++)
        {
            if (scheme1)
            {
                // Solve momentum equation on faces without the pressure gradient
                dhUdt = - hf*((F^Uf) & mesh.Sf())
                       - hf*magg*fvc::snGrad(h0)*mesh.magSf()
                       - fvc::flux(fvc::div(flux,U));
                flux = flux.oldTime()
                     + dt*((1-alpha)*dhUdt.oldTime() + alpha*dhUdt);

                // Create the pressure equation
                fvScalarMatrix hEqn
                (
                    fvm::ddt(h)
                  + fvc::div((1-alpha)*flux.oldTime() + alpha*flux)
                  - fvm::laplacian(sqr(alpha)*dt*magg*hf, h)
                );
                hEqn.solve();
                hf = fvc::interpolate(h);
                
                // Back substitutions
                if (alpha > 0) 
                {
                    flux += hEqn.flux()/alpha;
                    dhUdt += hEqn.flux()/(sqr(alpha)*dt);
                }

                // Update the velocity field on cell centres and faces
                U = fvc::reconstruct(flux/hf);
                Uf = fvc::interpolate(U);
            }         
            
            
            
            
            
            else if (scheme2) //1.2 with explicit mass advection
            {
                h.storePrevIter();
                hf.storePrevIter();
                U.storePrevIter();
                
                // Solve momentum equation on centres without the pressure gradient
                volVectorField hU = h.oldTime() * U.oldTime();
                surfaceScalarField phi_0
                (
                    fvc::flux(hU)
                );
                surfaceScalarField phi_1
                (
                    fvc::flux(h.prevIter() * U.prevIter())
                );
                              
                fvVectorMatrix huEqn
                (
                    fvm::ddt(h,U)
                  + (1-alpha)*fvc::div(phi_0,U.oldTime())
                  + alpha*fvm::div(phi_1,U)
                  + (1-alpha)*h.oldTime()*(Fc^U.oldTime())
                  + alpha*h.prevIter()*(Fc^U.prevIter())
                  + (1-alpha)*h.oldTime()*magg*fvc::grad(h.oldTime() + h0)
                  + alpha*h.prevIter()*magg*fvc::grad(h.prevIter() + h0)
                );
                huEqn.solve(); //to give u_prime
                
                
                hU = h.oldTime()*U.oldTime()
                   - dt*
                (
                    (1-alpha)*fvc::div(phi_0,U.oldTime())
                   + alpha*fvc::div(phi_1,U)
                   + (1-alpha)*h.oldTime()*(Fc^U.oldTime())
                   + alpha*h.prevIter()*(Fc^U.prevIter())
                   + (1-alpha)*h.oldTime()*magg*fvc::grad(h.oldTime()+h0)
                );
                
                U = hU/h.prevIter(); //to give u_prime_prime
                
                surfaceScalarField phi_prime
                (
                    fvc::flux(h.prevIter()*U)
                );
                
                // Create the pressure equation
                fvScalarMatrix hEqn
                (
                    fvm::ddt(h)
                  + fvc::div(phi_prime)
                  - fvm::laplacian(sqr(alpha)*dt*magg*hf.prevIter(), h)
                );
                hEqn.solve();
                hf = fvc::interpolate(h);
                
                // Back substitution
                if (alpha > 0) 
                {
                    U = U - alpha*dt*magg*fvc::grad(h + h0);
                }

                // Update the velocity field on cell faces
                //U = hU/h;
                Uf = fvc::interpolate(U);
                
            }
            
            
            
            
            
            
            else if (scheme3) //1.3 with implicit mass advection
            {
                h.storePrevIter();
                hf.storePrevIter();
                U.storePrevIter();
                
                // Solve momentum equation on centres without the pressure gradient
                volVectorField hU = h.oldTime() * U.oldTime();
                
                surfaceScalarField phi
                (
                    fvc::flux(hU)
                );
                surfaceScalarField phi_1
                (
                    fvc::flux(h.prevIter() * U.prevIter())
                );
                            
                              
                fvVectorMatrix huEqn
                (
                    fvm::ddt(h,U)
                  + (1-alpha)*fvc::div(phi,U.oldTime())
                  + alpha*fvm::div(phi_1,U)
                  + (1-alpha)*h.oldTime()*(Fc^U.oldTime())
                  + alpha*h.prevIter()*(Fc^U.prevIter())
                );
                huEqn.solve();
                
                U_S_new= fvc::flux(U);
                
                
                // Create the pressure equation
                fvScalarMatrix hEqn
                (
                    fvm::ddt(h)
                  + (1-alpha)*fvc::div(phi)
                  + (alpha)*fvm::div(U_S_new, h)
                  - fvc::laplacian(alpha*(1-alpha)*dt*magg*hf.oldTime(), h.oldTime())
                  - fvc::laplacian(alpha*(1-alpha)*dt*magg*hf.oldTime(), h0)
                  - fvm::laplacian(sqr(alpha)*dt*magg*hf.prevIter(), h)
                  - fvc::laplacian(sqr(alpha)*dt*magg*hf.prevIter(), h0)
                );
                hEqn.solve();
                hf = fvc::interpolate(h);
                
                // Back substitution
                if (alpha > 0) 
                {
                    hU = h*U - (1-alpha)*dt*magg*h.oldTime()*fvc::grad(h.oldTime()+h0) - alpha*dt*magg*h*fvc::grad(h+h0);
                }

                // Update the velocity field
                U = hU/h;
                Uf = fvc::interpolate(U);
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
