/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    BL2D

// ************************************************************************* //
// ****************                               ************************** //
// ****************   Alternative Linearisation   ************************** //
// ****************  Hilary's Method with Tristan ************************** //
// ************************************************************************* //   


Description
    Solves the Monge-Ampere equation to move a mesh based on a monitor
    function defined on the original mesh by using the Adaptive Linearisation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "monitorFunction.H"
#include "faceToPointReconstruct.H"


using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:



int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // Create the mesh to be moved
    fvMesh rMesh
    (
        Foam::IOobject
        (
            "rMesh", runTime.timeName(), runTime,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        )
    );

    // Open control dictionary
    IOdictionary controlDict
    (
        IOobject
        (
            args.executable() + "Dict", runTime.system(), runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    // The monitor funciton
    autoPtr<monitorFunction> monitorFunc(monitorFunction::New(controlDict));
    
    // The ratio of the finite different to the geometric Hessian
    const dimensionedScalar Gamma
    (
        controlDict.lookup("Gamma")
    );

    scalar conv = readScalar(controlDict.lookup("conv"));
    const int nCorr = readLabel(mesh.solutionDict().lookup("nCorrectors"));
       
    #include "createFields.H"


    Info << "Iteration = " << runTime.timeName()
         << " PABe = " << PABe.value() << endl;

    // Use time-steps instead of iterations to solve the Monge-Ampere eqn
    bool converged = false;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << flush << nl;

        phiBarLaplacian = fvc::laplacian(Phi);
        // Calculate the matrix: matrixA = 1+fvc::laplacian(phiBar)-Hessian
        tensor localA (0,0,0,0,0,0,0,0,0);
        scalar localAew = 0.0;
        bool print=true;
        forAll(matrixA, cellI)
        {
            
            matrixA[cellI] = diagTensor::one*(1+phiBarLaplacian[cellI]) - Hessian[cellI];
            matrixA[cellI].yy() = 1.0;
            //            localA = 0.5*(matrixA[cellI]+matrixA[cellI].T());
            localA = matrixA[cellI];
            localAew = eigenValues(localA)[0];
            if(localAew <= 0)
            {
                if(print)
                {
                    Info << "Minimum eigenvalue = " << localAew << endl;
                    print = false;
                }
                matrixA[cellI] = localA + (1.0e-5 - localAew)*diagTensor::one;
            }
            else
            {
                matrixA[cellI] = localA;
            }
        }
        // Calculate the source terms for the MA equation
        source = detHess - equiDistMean/monitorNew;

        // Setup and solve the MA equation to find Psi(t+1) 
        for (int iCorr = 0; iCorr < nCorr; iCorr++)
        {
            fvScalarMatrix PhiEqn
            (
                Gamma*fvm::laplacian(matrixA, Phi)
              - Gamma*fvc::laplacian(matrixA, Phi.oldTime())
              + source
            );
            PhiEqn.setReference(0, scalar(0));

            PhiEqn.solve();
        }

        // Calculate the gradient of phiBar at cell centres and on faces
        gradPhi = fvc::reconstruct(fvc::snGrad(Phi)*mesh.magSf());
        gradPhi.boundaryField()
            == (static_cast<volVectorField>(fvc::grad(Phi))).boundaryField();

        // Interpolate gradPhi (gradient of Phi) onto faces and correct the normal component
        gradPhif = fvc::interpolate(gradPhi);
        gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
                    *mesh.Sf()/mesh.magSf();

        // Map gradPhi onto vertices in order to create the new mesh
        pointVectorField gradPhiP
             = fvc::faceToPointReconstruct(fvc::snGrad(Phi));
        rMesh.movePoints(mesh.points() + gradPhiP);

        // finite difference Hessian of phiBar and its determinant
        Hessian = fvc::grad(gradPhif);
        forAll(detHess, cellI)
        {
            detHess[cellI] = det(diagTensor::one + Hessian[cellI]);
        }

       
        // Geometric version of the Hessian
        // detHess.internalField() =rMesh.V()/mesh.V();
        
        
        // map to or calculate the monitor function on the new mesh
        monitorR = monitorFunc().map(rMesh, monitor);
        monitorNew.internalField() = monitorR.internalField();
        monitorNew.correctBoundaryConditions();

        // The Equidistribution
        equiDist = monitorR*detHess;

        // mean equidistribution, c
        equiDistMean = fvc::domainIntegrate(detHess)/fvc::domainIntegrate(1/monitorNew);

       // The global equidistribution
       PABem = sum(equiDist)/mesh.nCells();
       PABe = pow((sum(pow((equiDist-PABem),2))/mesh.nCells()),0.5)/PABem;
       converged = PABe.value() < conv;

        Info << "Iteration = " << runTime.timeName()
             << " PABe = " << PABe.value() << endl;

        if (converged)
        {
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << nl <<endl;
            runTime.writeAndEnd();
        }
        runTime.write();
    }
    
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
