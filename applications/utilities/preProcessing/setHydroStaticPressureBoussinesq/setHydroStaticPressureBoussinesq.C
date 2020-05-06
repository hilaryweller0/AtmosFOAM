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
    setHydroStaticPressureBoussinesq

Description
    Find discretely hydrostatically balanced dynamic pressure P in balance with 
    a given buoyancy field b

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef CREATE_TIME
    #define CREATE_TIME createTime.H
#endif

#ifndef CREATE_MESH
    #define CREATE_MESH createMesh.H
#endif

#ifndef CREATE_FIELDS
    #define CREATE_FIELDS createFields.H
#endif

#ifndef CREATE_CONTROL
    #define CREATE_CONTROL createControl.H
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INCLUDE_FILE(X) INCLUDE_FILE2(X)
#define INCLUDE_FILE2(X) #X

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include INCLUDE_FILE(CREATE_TIME)
    #include INCLUDE_FILE(CREATE_MESH)
    
    #ifndef NO_CONTROL
    #include INCLUDE_FILE(CREATE_CONTROL)
    #endif
    
    const dictionary& itsDict = mesh.solutionDict().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
    //const scalar Ptol = itsDict.lookupOrDefault<scalar>("initPtol", SMALL);
    label pRefCell = mesh.solutionDict().lookupOrDefault<label>("pRefCell", 0);
    
    // Read in pressure and buoyancy
    Info<< "Reading field P\n" << endl;
    volScalarField P
    (
        IOobject
        (
            "P",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading in buoyancy, b\n" << endl;
    volScalarField b
    (
        IOobject
        (
            "b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // Surface scalar fields for hydrostatic boundary conditions
    const surfaceScalarField gradPcoeff
    (
        IOobject("gradPcoeff", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("1", dimless, scalar(1))
    );
    surfaceScalarField bf
    (
        "bf", 
        fvc::interpolate(b, "b")*mesh.Sf().component(2)
    );
    
    // Calculate and write out pressure that is in hydrostatic balance
    bool converged = false;
    
    Info << "Setting hydrostatic pressure " << endl;
    
    for (int iter = 0; iter <= maxIters && !converged; iter++)
    {
        Info << "Iteration: " << iter << endl;
        fvScalarMatrix PEqn
        (
            fvc::div(bf) - fvm::laplacian(P)
        );
        if (pRefCell >= 0) PEqn.setReference(pRefCell,0);
        //PEqn.solve();
        converged = PEqn.solve(P.name()).nIterations() == 0;
    }
    P.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

