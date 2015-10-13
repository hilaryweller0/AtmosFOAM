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
    setAdvectedTheta

Description
    Set theta based on an array of Brunt Vaisalla frequencies for different
    layers, transformed by the BTF function.  This creates an analytic solution
    of advecting a theta profile created with setTheta in a BTF velocity field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ExnerTheta.H"
#include "ThermalProfile.H"
#include "BTF.H"
#include "Mountain.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"

    ThermalProfile profile(envProperties, g, T0);

    IOdictionary dict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const SchaerExpMountain mountain(dict);
    const BTF btf(mountain, dict);
        
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
        const point p = mesh.C()[cellI];
        theta[cellI] = profile.thetaAt(btf.transform(p));
    }

    theta.write();

    Info<< "End\n" << endl;

    return 0;
}
