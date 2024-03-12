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
    timeMean

Description
    Calculates running time mean and outputs every input file

\*---------------------------------------------------------------------------*/


#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

template<class Type, template<class> class PatchField, class GeoMesh>
void timeAverage
(
    const instantList& timeDirs,
    const fileName& fieldName,
    const fileName& outFieldName,
    Time& runTime,
    const fvMesh& mesh
)
{
    GeometricField<Type,PatchField,GeoMesh> f
    (
        IOobject(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    GeometricField<Type,PatchField,GeoMesh> fMean(outFieldName, f);
    fMean.write();
    
    scalar startTime = timeDirs[0].value();
    scalar prevTime = startTime;
    scalar nowTime = startTime;
    scalar nextTime = startTime;
    
    for(label timeI = 1; timeI < timeDirs.size()-1; timeI++)
    {
        Info << "Averaging over time " << timeDirs[timeI].value() << endl;
        runTime.setTime(timeDirs[timeI], timeI);
        nowTime = timeDirs[timeI].value();
        nextTime = timeDirs[timeI+1].value();
        f = GeometricField<Type,PatchField,GeoMesh>
        (
            IOobject(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        fMean = 
        (
            (0.5*(nowTime+prevTime) - startTime)*fMean
           + 0.5*(nextTime - prevTime)*f
        )/(0.5*(nextTime+nowTime)-startTime);
        
        fMean.write();
        prevTime = nowTime;
    }
    if (timeDirs.size() > 1)
    {
        Info << "Averaging over time " << timeDirs.last().value() << endl;
        runTime.setTime(timeDirs.last(), timeDirs.size()-1);
        nowTime = timeDirs.last().value();
        f = GeometricField<Type,PatchField,GeoMesh>
        (
            IOobject(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        fMean = 
        (
            (0.5*(nowTime+prevTime) - startTime)*fMean
           + 0.5*(nowTime - prevTime)*f
        )/(nowTime-startTime);
    
        fMean.write();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("field");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    const fileName fieldName(args.args()[1].c_str());
    const fileName outFieldName=fieldName+"TimeMean";
    #include "createMesh.H"

    // Read in field for the first time step and initialise the average
    runTime.setTime(timeDirs[0], 0);
    
    IOobject fieldHeader(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ);

    if ( !fieldHeader.headerOk())
    {
        FatalErrorIn("timeMean") << "Cannot read " << fieldName
            << " from time directory " << runTime.timeName()
            << exit(FatalError);
    }
    else if (fieldHeader.headerClassName() == "volScalarField")
    {
        timeAverage<scalar, fvPatchField, volMesh>
        (
            timeDirs, fieldName, outFieldName, runTime, mesh
        );
    }
    else if (fieldHeader.headerClassName() == "surfaceScalarField")
    {
        timeAverage<scalar, fvsPatchField, surfaceMesh>
        (
            timeDirs, fieldName, outFieldName, runTime, mesh
        );
    }
    else if (fieldHeader.headerClassName() == "volVectorField")
    {
        timeAverage<vector, fvPatchField, volMesh>
        (
            timeDirs, fieldName, outFieldName, runTime, mesh
        );
    }
    else if (fieldHeader.headerClassName() == "surfaceVectorField")
    {
        timeAverage<vector, fvsPatchField, surfaceMesh>
        (
            timeDirs, fieldName, outFieldName, runTime, mesh
        );
    }
    else
    {
        FatalErrorIn("timeMean") << fieldName << " is of type "
            << fieldHeader.headerClassName()
            << " which is not supported by timeMean\n"
            << exit(FatalError);
    }
    
    Info << endl;
    return(0);
}


// ************************************************************************* //
