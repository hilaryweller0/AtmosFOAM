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
    writeuvw

Description
    Calculates the components of the vectorField specified

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("field");

#   include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    const word fieldName(args.additionalArgs()[0]);

    runTime.setTime(Times[startTime], startTime);

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        Info<< "    Reading vector " << fieldName << endl;

        IOobject fieldHeader
        (
            fieldName, runTime.timeName(), mesh, IOobject::MUST_READ
        );

        if ( !fieldHeader.headerOk())
        {
            FatalErrorIn("writeuvw") << "Cannot read " << fieldName
                << " from time directory " << runTime.timeName()
                << exit(FatalError);
        }
        else if (fieldHeader.headerClassName() == "volVectorField")
        {
            volVectorField U(fieldHeader, mesh);

            volScalarField U1(fieldName+"x", U.component(0));
            Info << "Writing " << fieldName+"x" << endl;
            U1.write();

            U1 == U.component(1);
            U1.rename(fieldName+"y");
            Info << "Writing " << fieldName+"y" << endl;
            U1.write();

            U1 == U.component(2);
            U1.rename(fieldName+"z");
            Info << "Writing " << fieldName+"z" << endl;
            U1.write();
        }
        else if (fieldHeader.headerClassName() == "surfaceVectorField")
        {
            surfaceVectorField U(fieldHeader, mesh);

            surfaceScalarField U1(fieldName+"x", U.component(0));
            Info << "Writing " << fieldName+"x" << endl;
            U1.write();

            U1 == U.component(1);
            U1.rename(fieldName+"y");
            Info << "Writing " << fieldName+"y" << endl;
            U1.write();

            U1 == U.component(2);
            U1.rename(fieldName+"z");
            Info << "Writing " << fieldName+"z" << endl;
            U1.write();
        }
        else
        {
            FatalErrorIn("writeuvw") << "Field type "
                << fieldHeader.headerClassName() << " of " << fieldName
                << " not supported" << exit(FatalError);
        }
    }

    return(0);
}


// ************************************************************************* //
