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

#include "fvCFD.H"
#include "argList.H"
#include "OFstream.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
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
    instantList timeDirs = timeSelector::select0(runTime, args);

    const word fieldName(args.additionalArgs()[0]);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    Foam::fvMesh mesh
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
    os << "#time mag RMS inf sum variance min max" << endl;
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
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
            
            scalar Vtot = 0;
            scalar l1 = 0;
            scalar l2 = 0;
            scalar li = 0;
            scalar l0 = 0;
            scalar min = 0.0;
            scalar max = 0.0;
            
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
                if (fi > max) max = fi;
                if (fi < min) min = fi;
            }
            l1 /= Vtot;
            l2 = Foam::sqrt(l2/Vtot);
            l0 /= Vtot;

            scalar variance = 0;
            forAll(sumCells, i)
            {
                label cellI = sumCells[i];
                scalar fi = f[cellI];
                scalar Vi = mesh.V()[cellI];
                variance += pow(fi - l0, 2) * Vi;
            }

            variance /= Vtot;
            
            os << runTime.timeName() << ' ' << l1 << ' ' << l2 << ' ' << li
               << ' ' << l0 << ' ' << variance << ' ' << min << ' ' << max << endl;
        }
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
