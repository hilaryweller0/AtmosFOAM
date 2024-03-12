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
    testDiv

Description
    Calculates div errors

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "linear.H"
#include "surfaceInterpolate.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensionedScalar L(dimLength, scalar(1));
    dimensionedScalar t(dimTime, scalar(1));

    Info << "Reading U\n" << endl;
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    surfaceScalarField phi("phi", t*linearInterpolate(U) & mesh.Sf());


    Info << "Reading T\n" << endl;
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    T = pow(mesh.C().component(0)/L, 3);

    volScalarField divT
    (
        "divT", fvc::div(phi,T)
    );
    
    volScalarField divExact = divT;
    divExact = 3*pow(mesh.C().component(0)/L, 2);
    
    volScalarField divError("divError", divT - divExact);
    divError.write();
    
    Info << "T = " << T.internalField();
    Info << "Interpolant = " << fvc::interpolate(T);
    Info << "divT = " << divT.internalField() << endl;
    Info << "divExact = " << divExact.internalField();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
