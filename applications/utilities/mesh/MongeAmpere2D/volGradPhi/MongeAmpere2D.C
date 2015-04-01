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
//#include "meshToPoint.H"
//#include "polyFit.H"
#include "meshToMesh.H"
#include "volPointInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    fvMesh rMesh
    (
        Foam::IOobject
        (
            "rMesh", runTime.timeName(), runTime,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        )
    );

    #include "createFields.H"
    
    // Use time-steps instead of iterations to solve the Monge-Ampere eqn
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << " "<<nl << flush;
        
        dimensionedScalar equiDistMean = 
             fvc::domainIntegrate(detHess-del2Phi)
                /fvc::domainIntegrate(1/monitorNew);

        fvScalarMatrix PhiEqn
        (
            fvm::laplacian(Phi)
          - equiDistMean/monitorNew
          - del2Phi
          + detHess
        );
        PhiEqn.setReference(0, scalar(0));
        PhiEqn.solve();

        gradPhi = fvc::grad(Phi);
        Hessian = fvc::grad(gradPhi);
        Hessian = 0.5*(Hessian + Hessian.T());
        del2Phi = fvc::laplacian(Phi);
        //del2Phi = fvc::div(fvc::grad(Phi));
        forAll(detHess, cellI)
        {
//            detHess[cellI] = det
//            (
//                diagTensor::one + Hessian[cellI]
//            );
            Hessian[cellI][4] = 1;
//            detHess[cellI] = 1 + del2Phi[cellI] + det(Hessian[cellI]);
            detHess[cellI] = 1 + tr(Hessian[cellI])-1 + det(Hessian[cellI]);
        }
        
        label cellI = 74;
        Info << "gradGrad = " << Hessian[cellI] << nl
             << "detHess = " << detHess[cellI] << nl
             << "det(H) = " << det(Hessian[cellI]) << nl
             << "tr(H) = " << tr(Hessian[cellI])-1 << nl
             << "del2Phi = " << del2Phi[cellI] << nl
             << "1+tr(H)+det(H) = " << tr(Hessian[cellI])+det(Hessian[cellI]) << nl
             << "1+del2+det(H) = " << 1 + del2Phi[cellI]+det(Hessian[cellI])
             << endl;
        
        //FatalErrorIn("") << exit(FatalError);
        
        // Geometric value of the Hessian
        //detHess.internalField() = rMesh.V()/mesh.V();
        // Correct the linear terms of detHess
        //detHess += del2Phi - tr(Hessian);
        
        // create the new mesh
        volPointInterpolation vpi(mesh);
        pointField gradPhiP = vpi.interpolate(gradPhi);
        rMesh.movePoints(mesh.points() + gradPhiP);
        //rMesh.write();

//        // Calculate the monitor function at the new locations
//        forAll(monitorNew, cellI)
//        {
//        // the new point of the cell
//            const point newPoint = mesh.C()[cellI] + gradPhi[cellI];
//            
//            // interpolate the monitor function onto the new point
//            meshToPoint<polyFit<1> > m2p(newPoint, mesh);
//            monitorNew[cellI] = m2p.interpolate(monitor);
//        }

        // map the monitor function onto the new mesh
        meshToMesh meshMap(mesh, rMesh, meshToMesh::imCellVolumeWeight, false);
        monitorR.internalField() = meshMap.mapSrcToTgt(monitor.internalField());
        monitorNew.internalField() = monitorR.internalField();
        
        equidist = monitorNew*detHess;
        
        Info << "Time = " << runTime.timeName()
             << " determinant goes from " << min(detHess).value()
             << " to " << max(detHess).value()
             << " equidist goes from "
             << min(equidist.internalField())
             << " to " << max(equidist.internalField())
             << " cell volumes go from " << min(rMesh.V()).value() << " to "
             << max(rMesh.V()).value() << " ratio = "
             << (max(rMesh.V())/min(rMesh.V())).value() << endl;

        if (min(detHess).value() < 0)
        {
            runTime.writeAndEnd();
        }
        runTime.write();
    }
    
    // write out new cell volumes
    volScalarField V
    (
        IOobject("V", runTime.timeName(), rMesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        rMesh, dimVol, rMesh.V(), monitor.boundaryField()
    );
    V.write();
    
//    Phi.write();
//    gradPhi.write();
//    monitorNew.write();
//    equidist.write();
//    detHess.write();
    Hessian.write();
//    del2Phi.write();
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
