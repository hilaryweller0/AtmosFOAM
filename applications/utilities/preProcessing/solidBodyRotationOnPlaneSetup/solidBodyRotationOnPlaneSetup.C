/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    solidBodyRotationOnPlaneSetup.C

Description
    Set the initial velocity field on faces, Uf, flux, phi and set an inital
    tracer field,T, for the solid body tracer advection test case described by
    Leonard, Lock and MacVean, MWR, 1996

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading initial conditions\n" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "solidBodyRotationOnPlaneSetupDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Constants in the specification of the velocity field
    const dimensionedVector rotationCentre(initDict.lookup("rotationCentre"));
    const dimensionedScalar spinSpeed(initDict.lookup("spinSpeed"));
    // Constanst in the specification of the tracer field
    const dimensionedScalar initTracerDist(initDict.lookup("initTracerDist"));
    const dimensionedScalar tracerInitAngle
    (
        dimensionedScalar(initDict.lookup("tracerInitAngle"))*pi/180
      + 2*runTime.time()*spinSpeed
    );
    const dimensionedScalar tracerRadius(initDict.lookup("tracerRadius"));
    
    // The initial tracer centre
    const dimensionedVector initialTracerCentre
         = rotationCentre + initTracerDist*vector
         (
            Foam::cos(tracerInitAngle.value()),
            Foam::sin(tracerInitAngle.value()),
            scalar(0)
         );
    
    // X and Y locations of the cell centres
    volScalarField X = mesh.C().component(vector::X);
    volScalarField Y = mesh.C().component(vector::Y);
    
    // Create initial flux and tracer fields
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    T = Foam::exp
    (
        -0.5*(magSqr(mesh.C() - initialTracerCentre)/sqr(tracerRadius))
    );
    T.write();
    
    // Circulate around all of the faces to set the flux from the 
    // streamfunction using Stokes theorem
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        point p0 = mesh.points()[f.last()];
        point p1 = mesh.points()[f.first()];
        point pmid = 0.5*(p0 + p1);
        vector streamFunction = spinSpeed.value()*vector(0,0,1)
                                *magSqr(pmid - rotationCentre.value());
        phi[faceI] = streamFunction & (p0 - p1);
        for(label ip = 1; ip < f.size(); ip++)
        {
            p0 = p1;
            p1 = mesh.points()[f[ip]];
            point pmid = 0.5*(p0 + p1);
            vector streamFunction = spinSpeed.value()*vector(0,0,1)
                                    *magSqr(pmid - rotationCentre.value());
            phi[faceI] += streamFunction & (p0 - p1);
        }
    }
    phi.write();
    
    // Set the velocity field from the flux
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        fvc::reconstruct(phi)
    );
    U.write();
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(U)
    );
    Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    Uf.write();
    
    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
