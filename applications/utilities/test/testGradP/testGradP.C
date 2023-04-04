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
    testGradP

Description
    Tests consistency of grad and snGrad

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "physicalProperties.H"
#include "fundamentalConstants.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "constTransport.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    surfaceScalarField snGradP
    (
        "snGradP",
        fvc::snGrad(P)
    );
    snGradP.write();

    volVectorField gradP
    (
        "gradP",
        fvc::reconstruct(snGradP*mesh.magSf()) - gradPexact
    );
    gradP.write();

    volVectorField gradP2("gradP2", fvc::grad(P) - gradPexact);
    gradP2.write();

    volVectorField gradP3("gradP3", fvc::grad(Pf) - gradPexact);
    gradP3.write();
    
    gradP2 += gradPexact;
    snGradP = fvc::interpolate(gradP2) & mesh.Sf()/mesh.magSf();
    volVectorField gradP4
    (
        "gradP4",
        fvc::reconstruct(snGradP*mesh.magSf()) - gradPexact
    );
    gradP4.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
