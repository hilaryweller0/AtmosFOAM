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
    makeHotBubble

Description
    Modify the initial (uniform) theta to include a warm bubble as defined by
    Bryan and Fritsch, MWR 2002

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    static scalar piby2 = 0.5*constant::mathematical::pi;

    const volVectorField& C = mesh.C();

    Info << "\nReading initialProperties" << endl;

    IOdictionary initialProperties
    (
        IOobject
        (
            "initialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar thetaPrime(initialProperties.lookup("thetaPrime"));
    vector bubbleCentre(initialProperties.lookup("bubbleCentre"));
    vector bubbleRadii(initialProperties.lookup("bubbleRadii"));

    Info<< "Reading theta_init\n" << endl;
    volScalarField theta_init
    (
        IOobject("theta_init", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info<< "Setting theta\n" << endl;
    volScalarField theta
    (
        IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
        theta_init
    );

    forAll(theta, celli)
    {
        scalar L = Foam::sqrt
        (
            sqr((C[celli].x() - bubbleCentre.x())/bubbleRadii.x())
          + sqr((C[celli].y() - bubbleCentre.y())/bubbleRadii.y())
          + sqr((C[celli].z() - bubbleCentre.z())/bubbleRadii.z())
        );

        if (L < 1)
        {
            theta[celli] += thetaPrime.value()*sqr(Foam::cos(piby2*L));
        }
    }

    Info<< "Writing theta\n" << endl;
    theta.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
