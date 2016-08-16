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
    setVelocityField

Description
    Sets the horizontal velocity field for the 
    Schar et al Mon. Wea. Rev., 130(10):2459-2480, 2002
    horizontal advection over orography test case

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "VelocityField.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
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

    autoPtr<VelocityProfile> velocityProfile(VelocityProfile::New(dict));
    const VelocityField velocityField(velocityProfile);

    Info << "Creating velocity field Uf" << endl;
    velocityField.applyTo(Uf);
    Uf.write();

    Info << "Creating flux field, phi" << endl;
    velocityField.applyTo(phi);
    phi.write();

    Info << "Calculating the divergence field to check that it is zero" << endl;
    volScalarField divu("divu", fvc::div(phi));
    divu.write();

    Info << "Correcting Uf based on the flux" << endl;
    Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    Uf.write();

    Info << "Creating velocity field U" << endl;
    velocityField.applyTo(U);
    U.write();
}

