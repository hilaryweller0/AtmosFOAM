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
    writes data in a x y z data table for the specified time step

Description
    

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

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
        << runTime.name() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.name(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    const word fieldName = args.args()[1];

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.name() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED
        );

        // initialise output file
        fileName outFile = args.rootPath() / args.caseName() / runTime.name();
        if (meshRegion != fvMesh::defaultRegion) outFile = outFile / meshRegion;
        outFile = outFile / fieldName + ".xyz";
        Info << "Writing file " << outFile << endl;
        OFstream os(outFile);
        os << "#x   y     z    " << fieldName << endl;


        if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == "volScalarField"
        )
        {
            volScalarField vf(fieldHeader, mesh);

            for(label celli = 0; celli < mesh.nCells(); celli++)
            {
                os << mesh.C()[celli][0] << "  " << mesh.C()[celli][1] << "  "
                   << mesh.C()[celli][2] << "  " << vf[celli] << endl;
            }
            if (!args.optionFound("noBoundary")) forAll(mesh.boundary(), patchI)
            {
                // Output the patch data if the patch has few faces than the 
                // mesh (ie if it covers new ground, ie not part of a 1D
                // or 2D domain)
                const fvPatch& patch = mesh.boundary()[patchI];
                const fvPatchField<scalar>& boundaryField = vf.boundaryField()[patchI];
                if(patch.size() < mesh.nCells() && !patch.coupled())
                    forAll(patch, faceI)
                {
                    os << patch.Cf()[faceI][0] << "  " << patch.Cf()[faceI][1] <<"  "
                       << patch.Cf()[faceI][2] << "  "
                       << boundaryField[faceI] << endl;
                }
            }
        }
        else if (fieldHeader.headerClassName() == "surfaceScalarField")
        {
            surfaceScalarField sf(fieldHeader, mesh);

            for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                os << mesh.Cf()[faceI][0] << "  " << mesh.Cf()[faceI][1] << "  "
                   << mesh.Cf()[faceI][2] << "  " << sf[faceI] << endl;
            }
        }
        else if (fieldHeader.headerClassName() == "surfaceVectorField")
        {
            surfaceVectorField vf(fieldHeader, mesh);

            for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                os << mesh.Cf()[faceI][0] << "  " << mesh.Cf()[faceI][1] <<"  "
                   << mesh.Cf()[faceI][2] << "  "
                   << vf[faceI][0] << "  "
                   << vf[faceI][1] << "  "
                   << vf[faceI][2] << endl;
            }
            forAll(mesh.boundary(), patchI)
            {
                const fvPatch& patch = mesh.boundary()[patchI];
                const fvsPatchField<vector>& boundaryField = vf.boundaryField()[patchI];
                forAll(patch, faceI)
                {
                    os << patch.Cf()[faceI][0] << "  " << patch.Cf()[faceI][1] <<"  "
                       << patch.Cf()[faceI][2] << "  "
                       << boundaryField[faceI][0] << "  "
                       << boundaryField[faceI][1] << "  "
                       << boundaryField[faceI][2] << endl;
                }
            }
        }
        else if (fieldHeader.headerClassName() == "volVectorField")
        {
            volVectorField vf(fieldHeader, mesh);

            for(label cellI = 0; cellI < mesh.nCells(); cellI++)
            {
                os << mesh.C()[cellI][0] << "  " << mesh.C()[cellI][1] << "  "
                   << mesh.C()[cellI][2] << "  "
                   << vf[cellI][0] << "  "
                   << vf[cellI][1] << "  "
                   << vf[cellI][2] << endl;
            }
            
            if (!args.optionFound("noBoundary"))  forAll(mesh.boundary(), patchI)
            {
                // Output the patch data if the patch has few faces than the 
                // mesh (ie if it covers new ground, ie not part of a 1D
                // or 2D domain)
                const fvPatch& patch = mesh.boundary()[patchI];
                const fvPatchField<vector>& boundaryField = vf.boundaryField()[patchI];
                if(patch.size() < mesh.nCells() && !patch.coupled())
                     forAll(patch, faceI)
                {
                    os << patch.Cf()[faceI][0] << "  " << patch.Cf()[faceI][1] <<"  "
                       << patch.Cf()[faceI][2] << "  "
                       << boundaryField[faceI][0] << "  "
                       << boundaryField[faceI][1] << "  "
                       << boundaryField[faceI][2] << endl;
                }
            }
        }
        else if (fieldHeader.headerClassName() == "volSymmTensorField")
        {
            volSymmTensorField tf(fieldHeader, mesh);
            
            for(label cellI = 0; cellI < mesh.nCells(); cellI++)
            {
                os << mesh.C()[cellI][0] << "  " << mesh.C()[cellI][1] << "  "
                   << mesh.C()[cellI][2] << "  "
                   << tf[cellI][0] << "  "
                   << tf[cellI][1] << "  "
                   << tf[cellI][2] << "  "
                   << tf[cellI][3] << "  "
                   << tf[cellI][4] << "  "
                   << tf[cellI][5] << endl;
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
