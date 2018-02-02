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
    setBalancedHeight

Description
    Find height field for the shallow water equations in balance with a
    given velocity field

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
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
    HodgeOps H(mesh);
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"

    const dictionary& itsDict = mesh.solutionDict().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Keep the domain averaged h fixed
    const dimensionedScalar hMean = fvc::domainIntegrate(h)
        /dimensionedScalar("", dimVol, gSum(mesh.V()));

    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Iteration " << iter << endl;
        hf = fvc::interpolate(h);
        V = hf * (Uf & H.delta());
        U = H.ddirToFlux(V);
        
        V = -dt*(H.delta() & fvc::interpolate(fvc::div(U,u), "Uf"));
        U = H.ddirToFlux(V - dt*hf*(two_dxOmega & Uf))
          - H.ddirToFluxCorr(dt*g*hf*fvc::snGrad(h)*H.magd());

        fvScalarMatrix hEqn
        (
            fvc::div(U) - fvm::laplacian(dt*g*hf, h, "h")
        );
        hEqn.setReference(0, h[0]);
        converged = hEqn.solve(mesh.solver(h.name())).nIterations() == 0;
        
        // Ensure the domain contains the correct mean h
        dimensionedScalar hMeanTmp = fvc::domainIntegrate(h)
            /dimensionedScalar("", dimVol, gSum(mesh.V()));
        h += hMean - hMeanTmp;
        Info << "h goes from " << min(h).value() << " to " << max(h).value()
             << " mean = " << fvc::domainIntegrate(h)
                              /dimensionedScalar("", dimVol, gSum(mesh.V()))
             << endl;
    }

    h.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

