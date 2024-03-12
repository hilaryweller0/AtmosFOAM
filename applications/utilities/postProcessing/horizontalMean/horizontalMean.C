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
    horizontalMean

Description
    Calculates the horizontal mean, standard deviation, min and max for
    a list of variables given in dictionary system/horizontalMeanDict.
    Also in the dictionary is a list of vertical level boundaries and for each
    field, need to say the density variable used to calculate volume average.
    Can use a subset of cells from a cellSet

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void mapFields
(
    const GeometricField<Type, fvPatchField, volMesh>& vfIn,
    GeometricField<Type, fvPatchField, volMesh>& vfOut,
    const label nConsec
)
{
    Info << "Mapping " << vfOut.name() << endl;

    label cellI = 0;
    forAll(vfOut, cellOut)
    {
        for(label i = 0; i < nConsec; i++)
        {
            vfOut[cellOut] += vfIn[cellI]*vfIn.mesh().V()[cellI];
            cellI++;
        }
        vfOut[cellOut] /= vfOut.mesh().V()[cellOut];
    }
    vfOut.correctBoundaryConditions();
    vfOut.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("targetCase");
    argList::validArgs.append("nConsecCells");
    argList::validArgs.append("fields");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    // Target case
    fileName targetPath = args[1];
    
    // Number of consecutive cells to average
    const label nConsec = readLabel(IStringStream(args[2])());

    // Initialize the set of selected fields from the command-line options
    wordList fields;
    IStringStream(args[3])() >> fields;

    // Target time
    Time runTimeTarget
    (
        Time::controlDictName,
        targetPath.path().toAbsolute(),
        fileName(targetPath.name())
    );
    
    // Target mesh
    fvMesh meshTarget
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    // Check that meshes are compatible
    label cellI = 0;
    forAll(meshTarget.V(), cellOut)
    {
        scalar volOut = 0;
        for(label i = 0; i < nConsec; i++)
        {
            volOut += mesh.V()[cellI];
            cellI++;
        }
        if (abs(volOut - meshTarget.V()[cellOut]) > SMALL)
        {
            FatalErrorIn("horizontalMean")
                 << " meshes are not compatible with nCconsec = " << nConsec
                 << exit(FatalError);
        }
    }

    // map fields for all times
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        runTimeTarget.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        forAll(fields, fieldI)
        {
            IOobject fieldHeader(fields[fieldI], runTime.timeName(), mesh,
                                IOobject::MUST_READ);
        
            if ( !fieldHeader.headerOk())
            {
                FatalErrorIn("horizontalMean") << "Cannot read "
                    << fields[fieldI] << " from time directory "
                    << runTime.timeName() << exit(FatalError);
            }
            else if (fieldHeader.headerClassName() == "volScalarField")
            {
                volScalarField vfIn(fieldHeader, mesh);
                volScalarField vfOut
                (
                    IOobject(fields[fieldI], runTimeTarget.timeName(), meshTarget),
                    meshTarget,
                    dimensionedScalar(fields[fieldI], vfIn.dimensions(), scalar(0))
                );
                
                mapFields<scalar>(vfIn, vfOut, nConsec);
            }
            else if (fieldHeader.headerClassName() == "volVectorField")
            {
                volVectorField vfIn(fieldHeader, mesh);
                volVectorField vfOut
                (
                    IOobject(fields[fieldI], runTimeTarget.timeName(), meshTarget),
                    meshTarget,
                    dimensionedVector(fields[fieldI], vfIn.dimensions(), vector::zero)
                );
                
                mapFields<vector>(vfIn, vfOut, nConsec);
            }
            else
            {
                FatalErrorIn("horizontalMean") << fields[fieldI] << " is of type "
                    << fieldHeader.headerClassName()
                    << " which is not handled by horizontalMean.\n"
            << "horizontalMean only handles volScalarFields and volVectorFields"
                    << exit(FatalError);
            }
        }
    }
    return(0);
}


// ************************************************************************* //
