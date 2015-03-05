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
    setHorizontalVelocityField

Description
    Sets the horizontal velocity field for the 
    Schar et al Mon. Wea. Rev., 130(10):2459-2480, 2002
    horizontal advection over orography test case

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "VelocityField.H"
#include "HorizontalVelocityProfile.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    Info << "Reading velocityFieldDict" << endl;

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

    const scalar u0(readScalar(dict.lookup("maxVelocity")));
    const scalar z1(readScalar(dict.lookup("zeroVelocityHeight")));
    const scalar z2(readScalar(dict.lookup("maxVelocityHeight")));

    const HorizontalVelocityProfile profile(u0, z1, z2);
    const VelocityField velocityField(profile);

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    Info << "Creating velocity field Uf" << endl;
    velocityField.applyTo(Uf);

    Uf.write();

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    Info << "Creating velocity field U" << endl;
    velocityField.applyTo(U);

    U.write();
}
