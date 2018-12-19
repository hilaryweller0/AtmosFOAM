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

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()
    #define rMesh mesh
    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
        )
    );

    // Read in and create variables and fields
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"

    const dictionary& itsDict = mesh.solutionDict().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Keep the domain averaged h fixed
    const dimensionedScalar meshVol("meshVol", dimVol, gSum(mesh.V()));
    const dimensionedScalar hMean = fvc::domainIntegrate(h)/meshVol;

    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Iteration " << iter << endl;
        hf = fvc::interpolate(h);
        phi = hf*(Uf & mesh.Sf());
        
        surfaceScalarField phiTmp = -dt*fvc::flux(fvc::div(phi,U))
                                    -dt*hf*fvc::flux(twoOmega ^ U);
        
        if (withMountain)
        {
            phiTmp -= dt*g*hf*fvc::snGrad(h0)*mesh.magSf();
        }
        
        fvScalarMatrix hEqn
        (
            fvc::div(phiTmp) - fvm::laplacian(dt*g*hf, h, "laplacian(h)")
        );
        hEqn.setReference(0, h[0]);
        converged = hEqn.solve(h.name()).nIterations() == 0;

        // Ensure the domain contains the correct mean h
        dimensionedScalar hMeanTmp = fvc::domainIntegrate(h)/meshVol;
        h += hMean - hMeanTmp;
        Info << "h goes from " << min(h).value() << " to " << max(h).value()
             << " mean = " << fvc::domainIntegrate(h)/meshVol
             << endl;
    }

    h.write();
    
    volVectorField hU("hU", h*U);
    hU.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

