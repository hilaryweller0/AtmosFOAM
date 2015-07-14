/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    perturbField

Description
    Sets up initial conditions for theta for the perturbation specified
    by skamarock-klemp1994

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    IOdictionary dict
    (
        IOobject
        (
            "perturbFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar d_theta_0(readScalar(dict.lookup("amplitude")));
    const scalar H(readScalar(dict.lookup("height")));
    const scalar x_c(readScalar(dict.lookup("x_centre")));
    const scalar a(readScalar(dict.lookup("halfWidth")));

    forAll(theta, cellI)
    {
        const point& p = mesh.C()[cellI];
        theta[cellI] += d_theta_0 * Foam::sin(constant::mathematical::pi * p.z() / H) / 
            (1 + pow(p.x() - x_c, 2) / (a*a));
    };
    theta.write();

    forAll(thetaf, faceI)
    {
        const point& p = mesh.Cf()[faceI];
        thetaf[faceI] += d_theta_0 * Foam::sin(constant::mathematical::pi * p.z() / H) / 
                (1 + pow(p.x() - x_c, 2) / (a*a));
    }
    thetaf.write();
}
