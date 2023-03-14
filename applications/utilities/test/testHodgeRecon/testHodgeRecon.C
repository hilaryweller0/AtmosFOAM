/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    testHodgeRecon

Description
    Test that the Hodge reconstruction operators can reconstruct uniform 
    velocity fields on distorted meshes

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
#   include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    Info << "Mesh has normal direction" << flush;
    vector meshNormal = 0.5*(Vector<label>(1,1,1)-mesh.geometricD());
    meshNormal -= 2*meshNormal[1]*vector(0.,1.,0.);
    Info << meshNormal << endl;

    // Declare the Hodge operator class for this mesh
    HodgeOps H(mesh);

    // Uniform velocity field to reconstruct
    surfaceVectorField Ufone
    (
        IOobject("Ufone", runTime.timeName(), mesh),
        //mesh,
        //dimensionedVector("Ufone", dimless, vector::one - meshNormal)
        mesh.Cf()
    );
    Ufone.write();
    
    surfaceVectorField Ufrecon
    (
        "Ufrecon",
        fvc::interpolate(H.reconstructd(Ufone & H.delta()))
        //fvc::interpolate(fvc::reconstruct(Ufone & mesh.Sf()))
    );
    // Ensure that Ufrecon is correct in the d direction
    Ufrecon += ((Ufone-Ufrecon) & H.delta())*H.delta()/sqr(H.magd());
    Ufrecon.write();
    
    surfaceVectorField UfDiff("UfDiff", Ufone - Ufrecon);
    UfDiff.write();
    
    surfaceScalarField magDiff("magUfDiff", mag(UfDiff));
    scalar Ufmax = max(mag(Ufone.internalField())).value();
    Info << "Error fraction goes from "
         << min(magDiff.internalField()).value()/Ufmax << " to "
         << max(magDiff.internalField()).value()/Ufmax << endl;

    // Error in reconstrucing the flux
    surfaceScalarField fluxRecon("fluxRecon", H.ddirToFlux(Ufone & H.delta()));
    scalar maxFlux = max(mag(Ufone& mesh.Sf())).value();
    surfaceScalarField fluxReconError
    (
        "fluxReconError", (fluxRecon - (Ufone& mesh.Sf()))/maxFlux
    );
    fluxReconError.write();
    Info << "Flux error fraction goes from "
         << min(fluxReconError.internalField()).value() << " to "
         << max(fluxReconError.internalField()).value() << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
