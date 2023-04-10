/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    writes components of wind in local Cartesian coords on the sphere

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "polarPoint.H"
#include "sphericalVector.H"
#include "sphericalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    Foam::argList::validArgs.append("fieldName");
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    Foam::argList::addOption("noBoundary");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    
    const word fieldName = args.args()[1];
    
    // Read earth radius from system/extrudeMeshDict
    IOdictionary dict
    (
        IOobject("extrudeMeshDict", runTime.system(), mesh, IOobject::MUST_READ)
    );
    dictionary eDict = dict.subDict("linearRadialCoeffs");
    const scalar earthRadius(readScalar(eDict.lookup("Rsurface")));

    // Calculate spherical mesh data
    sphericalMeshData spherical(mesh, earthRadius);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED
        );

        // initialise output file
        fileName outFile = args.rootPath() / args.caseName() / runTime.timeName();
        if (meshRegion != fvMesh::defaultRegion) outFile = outFile / meshRegion;
        outFile = outFile / fieldName + ".latLon";
        Info << "Writing file " << outFile << endl;
        OFstream os(outFile);
        os << "#lon   lat     z    " << fieldName << endl;

        if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == "volVectorField"
        )
        {
            volVectorField vf(fieldHeader, mesh);

            for(label celli = 0; celli < mesh.nCells(); celli++)
            {
                point xyz = spherical.cellCentres()[celli];
                polarPoint C = convertToPolar(xyz, 360, earthRadius);
                sphericalVector wind = convertToLocal(xyz, vf[celli]);
                
                os << C.Lon() << "  " << C.Lat() << "  " << C.Z() << "  "
                   << wind.v[0] << "  " << wind.v[1] << "  " << wind.v[2] << endl;
            }
        }
        else if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == "surfaceVectorField"
        )
        {
            surfaceVectorField vf(fieldHeader, mesh);

            for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
            {
                point xyz = spherical.faceCentres()[facei];
                polarPoint Cf = convertToPolar(xyz, 360, earthRadius);
                sphericalVector wind = convertToLocal(xyz, vf[facei]);
                
                os << Cf.Lon() << "  " << Cf.Lat() << "  " << Cf.Z() << "  "
                   << wind.v[0] << "  " << wind.v[1] << "  " << wind.v[2] << endl;
            }
        }
        else
        {
            Info << "No " << fieldName << endl;
        }
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
