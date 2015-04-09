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
    MongeAmpere2D

Description
    Solves the Monge-Ampere equation to move a mesh based on a monitor
    function defined on the original mesh.

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
    // (HessianVolumeRatio=1 => geometric Hessian,
    //  HessianVolumeRatio=0 => finite difference Hessian)
    const scalar HessianVolumeRatio
    (
        readScalar(controlDict.lookup("HessianVolumeRatio"))
    );

    dimensionedScalar Vtot = sum(mesh.V());

    // boost the Laplacian to stabilise
    scalar boostLaplacian(readScalar(controlDict.lookup("boostLaplacian")));

    #include "createFields.H"

    // Use time-steps instead of iterations to solve the Monge-Ampere eqn
    bool converged = false;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << flush << nl;

        // Calculate the source terms for the MA equation
        source = detHess - equiDistMean/monitorNew;

        // calculate boostLaplacian to achieve stability
        boostLaplacian = max
        (
            boostLaplacian,
            4*max(scalar(0.25), max(mag(source.internalField())))
        );

        // Setup and solve the MA eqn (expressed as a Poisson equation)
        fvScalarMatrix PhiEqn
        (
            boostLaplacian*fvm::laplacian(Phi)
          - boostLaplacian*fvc::laplacian(Phi)
          + source
        );
        PhiEqn.setReference(0, scalar(0));
        solverPerformance sp = PhiEqn.solve();
        converged = sp.nIterations() <= 1;

        // Calculate the gradient of Phi at cell centres and on faces
        gradPhi = fvc::reconstruct(fvc::snGrad(Phi)*mesh.magSf());
        gradPhi.boundaryField()
            == (static_cast<volVectorField>(fvc::grad(Phi))).boundaryField();
        // Interpolate gradPhi onto faces and correct the normal component
        gradPhif = fvc::interpolate(gradPhi);
        gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
                    *mesh.Sf()/mesh.magSf();

        // Map gradPhi onto vertices in order to create the new mesh
        pointVectorField gradPhiP
             = fvc::faceToPointReconstruct(fvc::snGrad(Phi));
        rMesh.movePoints(mesh.points() + gradPhiP);

        // finite difference Hessian and its determinant
        Hessian = fvc::grad(gradPhif);
        forAll(detHess, cellI)
        {
            detHess[cellI] = det(diagTensor::one + Hessian[cellI]);
        }
        // Geometric version of the Hessian
        volRatio.internalField() =rMesh.V()/mesh.V();
        
        // combine the finite difference and geometric Hessian
        detHess = (1-HessianVolumeRatio)*detHess
                + HessianVolumeRatio*volRatio;

        // map to or calculate the monitor function on the new mesh
        monitorR = monitorFunc().map(rMesh, monitor);
        monitorNew.internalField() = monitorR.internalField();
        monitorNew.correctBoundaryConditions();

        // mean equidistribution
        equiDistMean = fvc::domainIntegrate(detHess)
                       /fvc::domainIntegrate(1/monitorNew);

        // Smooth the monitor function for the source term
        //monitorNew += 0.25*fvc::laplacian(sqr(1/mesh.deltaCoeffs()), monitorNew);
        //monitorNew += 0.25*fvc::laplacian(sqr(1/mesh.deltaCoeffs()), monitorNew);
        
        Info << "Time = " << runTime.timeName()
             << " source goes from "
             << min(source.internalField())
             << " to " << max(source.internalField())
             << " cell volumes go from " << min(rMesh.V()).value() << " to "
             << max(rMesh.V()).value()
             << " boostLaplacian = " << boostLaplacian << endl;

        if (converged)
        {
            runTime.writeAndEnd();
        }
        runTime.write();
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
