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
    exnerFoamH

Description
    Transient Solver for buoyant, inviscid, incompressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and V (flux in d direction)

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "atmosphere.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "orthogonalBoundaries.H"
//    #include "readEnvironmentalProperties.H"
//    #include "readThermoProperties.H"
//    #include "readThermoPropertiesMoist.H"
//    HodgeOps H(mesh);
//    surfaceScalarField gd("gd", g & H.delta());
//    #define dt runTime.deltaT()
    #include "createFields.H"

    runTime++;
    for(label ip = 0; ip < atmos.size(); ip++)
    {
        Info << "Species " << ip << endl;
        atmos[ip].gas().rho().write();
        atmos[ip].liquid().v().write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
