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
    Creates a computedThetaf field from an existing thetaf by transforming
    it via a buoyancy

Description
    

\*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include "fvCFD.H"
#include "OFstream.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"

    Info<< "Reading field thetaf\n" << endl;
    surfaceScalarField thetaf
    (
        IOobject
        (
            "thetaf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField bf
    (
        IOobject
        (
            "bf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thetaf * gUnitNormal
    );
    bf.write();

    volVectorField b = fvc::reconstruct(bf * mesh.magSf());

    surfaceScalarField computedThetaf("computedThetaf", mag(bf) + fvc::interpolate(b & ghat) * (1.0 - mag(gUnitNormal)));
    computedThetaf.write();

    volScalarField computedTheta
    (
        IOobject
        (
            "computedTheta",
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        b & ghat
    );
    computedTheta.write();

    return EXIT_SUCCESS;
}
