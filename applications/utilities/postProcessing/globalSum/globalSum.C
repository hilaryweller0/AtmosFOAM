/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    Werrors

Description
    Calculates the Williamson et al error measures based on differences between
    target and reference solutions. Can use a subset of cells from a cellSet

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "argList.H"
#include "OFstream.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("field");
    argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    argList::addOption
    (
        "cellSet", "cellSetName", "only calculate sums for a subset of cells"
    );

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    const word fieldName(args.additionalArgs()[0]);

    runTime.setTime(Times[startTime], startTime);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    Foam::fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // List of cells to sum
    labelList sumCells;
    if (args.optionFound("cellSet"))
    {
        const word cellSetName(args.optionRead<string>("cellSet"));
        cellSet cells(mesh, cellSetName);
        sumCells = cells.toc();
    }
    else
    {
        // Select all cells
        sumCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            sumCells[cellI] = cellI;
        }
    }

    // initialise diagnostics file
    fileName outFile = args.optionFound("region") ?
        args.rootPath() / args.caseName() /"globalSum"+meshRegion+fieldName+".dat"
      : args.rootPath() / args.caseName() /"globalSum"+fieldName+".dat";
    Info << "Writing global sums to " << outFile << endl;
    OFstream os(outFile);
    os << "#time mag RMS inf sum" << endl;
    
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject header
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (!header.headerOk())
        {
            Info << "No " << fieldName << endl;
        }
        else if(header.headerClassName() == "volScalarField")
        {
            const volScalarField f(header, mesh);

//            dimensionedScalar Vtot = sum(mesh.V());
//            scalar l1 = (sum(mag(f)*mesh.V())/Vtot).value();
//            scalar l2 = (Foam::sqrt(sum(sqr(f)*mesh.V())/Vtot)).value();
//            scalar li = (max(mag(f))).value();
//            scalar l0 = (sum(f*mesh.V())/Vtot).value();
            
            scalar Vtot = 0;
            scalar l1 = 0;
            scalar l2 = 0;
            scalar li = 0;
            scalar l0 = 0;
            
            forAll(sumCells, i)
            {
                label cellI = sumCells[i];
                scalar fi = f[cellI];
                scalar Vi = mesh.V()[cellI];
                Vtot += Vi;
                l1 += mag(fi)*Vi;
                l2 += sqr(fi)*Vi;
                if (mag(fi) > li) li = mag(fi);
                l0 += fi*Vi;
            }
            l1 /= Vtot;
            l2 = Foam::sqrt(l2/Vtot);
            l0 /= Vtot;
            
            os << runTime.timeName() << ' ' << l1 << ' ' << l2 << ' ' << li
               << ' ' << l0 << endl;
        }
//        else if(header.headerClassName() == "surfaceScalarField")
//        {
//            const surfaceScalarField f(header, mesh);
//            
//            scalar Vtot = sum(mesh.V()).value();
//            scalar l1 = sum(TRiSK::volumeIntegrate(mag(f)))/Vtot;
//            scalar l2 = Foam::sqrt(sum(TRiSK::volumeIntegrate(sqr(f)))/Vtot);
//            scalar li = (max(mag(f))).value();
//            scalar l0 = sum(TRiSK::volumeIntegrate(f))/Vtot;
//            
//            os << runTime.timeName() << ' ' << l1 << ' ' << l2 << ' ' << li
//               << ' ' << l0 << endl;
//        }
        else
        {
            FatalErrorIn("globalSum") << "Field type "
                << header.headerClassName() << " of " << fieldName
                << " not supported" << exit(FatalError);
        }
    }
    return(0);
}


// ************************************************************************* //
