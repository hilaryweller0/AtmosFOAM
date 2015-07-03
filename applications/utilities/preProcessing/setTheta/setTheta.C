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
    setTheta

Description
    Set theta based on an array of Brunt Vailsalla frequencies for different
    layers

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ExnerTheta.H"
#include "ThermalProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"

    ThermalProfile profile(envProperties, g, T0);
        
    Info<< "Reading theta_init\n" << endl;
    volScalarField theta_init
    (
        IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    // theta
    Info<< "Creating theta\n" << endl;
    volScalarField theta
    (
        IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
        theta_init
    );

    forAll(theta, cellI)
    {
        theta[cellI] = profile.thetaAt(mesh.C()[cellI]);
    }
    forAll(theta.boundaryField(), patchI)
    {
        fvPatchField<scalar>& thetap = theta.boundaryField()[patchI];
        forAll(thetap, facei)
        {
            thetap[facei] = profile.thetaAt(mesh.C().boundaryField()[patchI][facei]);
        }
    }
    theta.write();

    // thetaf
    surfaceScalarField thetaf("thetaf", fvc::interpolate(theta));
    forAll(thetaf, faceI)
    {
        thetaf[faceI] = profile.thetaAt(mesh.Cf()[faceI]);
    }
    forAll(thetaf.boundaryField(), patchI)
    {
        fvsPatchField<scalar>& thetap = thetaf.boundaryField()[patchI];
        forAll(thetap, faceI)
        {
            thetap[faceI] = profile.thetaAt(mesh.Cf().boundaryField()[patchI][faceI]);
        }
    }
    thetaf.write();

    volScalarField BruntFreq
    (
        IOobject("BruntFreq", runTime.timeName(), mesh),
        Foam::sqrt(-(g & fvc::grad(thetaf))/theta)
    );
    BruntFreq.write();

    surfaceScalarField BruntFreq2f
    (
        IOobject("BruntFreq2f", runTime.timeName(), mesh),
        -(g & mesh.delta())/mag(mesh.delta())*fvc::snGrad(theta)/thetaf
    );
    BruntFreq2f.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

