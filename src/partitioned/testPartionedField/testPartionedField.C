/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    testPartionedField

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    wordList partNames(2);
    partNames[0] = "stable";
    partNames[1] = "buoyant";

    Info << "Reading in partionedVolScalarFraction sigma" << endl;
    partitionedVolScalarFraction sigma
    (
        partNames, 
        IOobject("sigma", runTime.timeName(), mesh),
        mesh
    );
    
    Info << "Reading in partionedVolScalarField rho" << endl;
    partitionedVolScalarField rho
    (
        partNames, 
        IOobject("rho", runTime.timeName(), mesh),
        mesh,
        sigma
    );
    
    partitionedVolScalarFraction sigmaRho = rho.fraction();

    rho = sigmaRho.field(sigma);

    runTime++;
    
    sigma.write();
    rho.write();
    sigmaRho.write();
    
    //sigma = sigmaRho;
    sigma.write();

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
