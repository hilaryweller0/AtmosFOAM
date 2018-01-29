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

    Info << "Reading in partionedVolScalarField sigma" << endl;
    partitionedVolScalarField sigma
    (
        "sigma", partNames, mesh, runTime.timeName()
    );
    
    Info << "Reading in partionedVolScalarField rho" << endl;
    partitionedVolScalarField rho
    (
        "rho", partNames, mesh, runTime.timeName(), sigma
    );
    
    partitionedVolScalarField sigmaRho = rho.fraction();

    rho = sigmaRho.divideBy(sigma);

    runTime++;
    
    sigmaRho.storeTime();
    
    sigma.write();
    rho.write();
    sigmaRho.write();
    
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
