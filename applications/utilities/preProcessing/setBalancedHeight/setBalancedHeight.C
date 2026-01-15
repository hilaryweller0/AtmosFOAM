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
    setBalancedHeightRC

Description
    Find height field for the shallow water equations in balance with a
    given velocity field with velocity, U, defined at cell centres

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
#include "fvcVolumeIntegrate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "readEarthProperties.H"
    #include "createFields.H"

    const dictionary& itsDict = mesh.solution().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Keep the domain averaged h fixed
    const dimensionedScalar meshVol(dimVolume, gSum(mesh.V()));
    const dimensionedScalar hMean = fvc::domainIntegrate(h)/meshVol;

    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Iteration " << iter << endl;
        
        hf = fvc::interpolate(h);
        phi = hf*fvc::flux(U);
        dhUdt = -fvc::div(phi, U) -h*(F ^ U);
        phi = dt*(fvc::flux(dhUdt) - magg*hf*fvc::snGrad(h0)*mesh.magSf());
        
        fvScalarMatrix hEqn
        (
            fvc::div(phi) - fvm::laplacian(dt*magg*hf, h, "laplacian(h)")
        );
        hEqn.setReference(0, h[0]);
        converged = hEqn.solve(h.name()).nIterations() == 0;

        // Ensure the domain contains the correct mean h
        dimensionedScalar hMeanTmp = fvc::domainIntegrate(h)/meshVol;
        h += hMean - hMeanTmp;
        Info << "h goes from " << min(h).value() << " to " << max(h).value()
             << " mean = " << (fvc::domainIntegrate(h)/meshVol).value()
             << endl;
    }

    h.write();
    
    // How close is the initial gradient wind balance
    hf = fvc::interpolate(h);
    phi = hf*fvc::flux(U);
    ghGradh = magg*hf*fvc::snGrad(h+h0)*mesh.magSf();
    volScalarField gradWind("gradWind", mag(fvc::div(phi,U) + h*(F^U)));
    gradWind.write();
    volScalarField ghGradhC("ghGradh", mag(fvc::reconstruct(ghGradh)));
    ghGradhC.write();
    volScalarField imbalance("imbalance", gradWind-ghGradhC);
    imbalance.write();
    
    volVectorField hU("hU", h*U);
    hU.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

