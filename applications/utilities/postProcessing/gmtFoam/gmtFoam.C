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
    gmtFoam

Description
    Plots data on a 2d OpenFOAM case or on a patch of a 3d case using gmt

\*---------------------------------------------------------------------------*/

//#include "meshWithDual.H"
#include "fvCFD.H"
#include "unitVectors.H"
#include "pointFields.H"
#include "OFstream.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "mathematicalConstants.H"
#include "stringScalar.H"
#include "IOmanip.H"

#include "setRegion.H"
#include "FieldToPlot.H"

#include "systemVerbose.H"
#include "polarPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("gmtDict");

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

    runTime.setTime(Times[startTime], startTime);

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh/*WithDual*/ mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
    
#   include "readGmtDict.H"

    const label patchi = patchName == "" ? -1 :
                         mesh.boundaryMesh().findPatchID(patchName);

    const bool plotAllCells = (patchi == -1);

    const polyPatch& plotPatch = mesh.boundaryMesh()[max(patchi, 0)];

    const scalar radToDeg = 180./constant::mathematical::pi;

#   include "checkProjection.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        // create output file
        fileName epsFile = epsFileName != args.additionalArgs()[0]
                         ? epsFileName + ".ps"
                         : args.rootPath() / args.caseName() / runTime.timeName()
                           / epsFileName + ".ps";
        fileName pdfFile = epsFileName != args.additionalArgs()[0]
                         ? epsFileName + ".pdf"
                         : args.rootPath() / args.caseName() / runTime.timeName()
                           / epsFileName + ".pdf";

        Info << "Creating " << pdfFile << endl;

        if ( !overlay)
        {
            systemCall = "echo 0 0 | psxy -R" + region + " -J"
                              + projection + " -A -K > " + epsFile;
            systemVerbose(systemCall);
        }

        // Plot all the required fields
        forAll(fieldsToPlot, ifield)
        {
            IOobject fieldHeader
            (
                fieldsToPlot[ifield].name(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (!fieldHeader.headerOk())
            {
                fieldHeader = IOobject
                (
                    fieldsToPlot[ifield].name(),
                    runTime.findInstance(mesh.meshDir(), "points"),
                    mesh.meshSubDir,
                    mesh,
                    IOobject::MUST_READ
                );
            }
            if
            (
                !fieldHeader.headerOk() &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::MESH &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::MESHPOINTS &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::MESHCENTRES &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::ADVECTED_CONTOURS &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::NUMBERED &&
                fieldsToPlot[ifield].plotType() != FieldToPlot::WRITECONTOURS
            )
            {
                Info << "No " << fieldsToPlot[ifield].name() << endl;
            }
            else
            {
                Info << "Plotting " << fieldsToPlot[ifield].name()
                     << " of type " << fieldHeader.headerClassName()
                     << " using " << fieldsToPlot[ifield].plotType() << endl;
            
                scalar colourMin = fieldsToPlot[ifield].min();
                scalar colourMax = fieldsToPlot[ifield].max();
                scalar colourStep = fieldsToPlot[ifield].delta();

                if (fieldHeader.headerClassName() == "volScalarField")
                {
                    volScalarField vf(fieldHeader, mesh);

                    switch (fieldsToPlot[ifield].plotType())
                    {
                        case FieldToPlot::FILLED_CONTOURS:
                            if (plotAllCells)
                            {
#                               include "allFilledVolContours.H"
                            }
                            else
                            {
#                               include "filledVolContours.H"
                            }
                            break;
                        case FieldToPlot::SOLID_CONTOURS:
                            if (plotAllCells)
                            {
#                               include "allVolContours.H"
                            }
                            else
                            {
#                               include "volContours.H"
                            }
                            break;
                        case FieldToPlot::DASHED_CONTOURS:
                            if (plotAllCells)
                            {
#                               include "allVolContours.H"
                            }
                            else
                            {
#                               include "volContours.H"
                            }
                            break;
                        case FieldToPlot::RAW_VALUES:
                            if (plotAllCells)
                            {
#                               include "allRawVolValues.H"
                            }
                            else
                            {
#                               include "rawVolValues.H"
                            }
                            break;
                        case FieldToPlot::WRITECONTOURS:
                            {
#                           include "writeContours.H"
                            }
                            break;
                        default:
                            FatalErrorIn("gmtFoam")
                            << fieldsToPlot[ifield].plotType()
                            << " is not a valid plot type for volScalarField "
                            << fieldsToPlot[ifield].name()
<< ". Valid types are filledContours, solidContours, dashedContours, writeContours or rawValues"
                            << exit(FatalError);
                    }
                }
                else if (fieldHeader.headerClassName()=="pointScalarField")
                {
                    if
                    (
                fieldsToPlot[ifield].plotType() == FieldToPlot::FILLED_CONTOURS
                        && !plotAllCells
                    )
                    {
#                       include "filledPointContours.H"
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                        << fieldsToPlot[ifield].plotType()
                        << " is not a valid plot type for pointScalarField "
                        << fieldsToPlot[ifield].name()
                        << ". Valid types are filledContours and "
                        << "you must specify a valid patch to plot on"
                        << exit(FatalError);
                    }
                }
                else if (fieldHeader.headerClassName()=="volVectorField")
                {
                    volVectorField U(fieldHeader, mesh);

                    scalar vectorScale = fieldsToPlot[ifield].vectorScale();
                    label vectorFreq   = fieldsToPlot[ifield].vectorFreq();

                    switch (fieldsToPlot[ifield].plotType())
                    {
                        case FieldToPlot::VECTORS:
                            if (plotAllCells)
                            {
#                               include "allVolVectors.H"
                            }
                            else
                            {
#                               include "volVectors.H"
                            }
                            break;
                        case FieldToPlot::VECTOR_END_POINTS:
                            if (plotAllCells)
                            {
#                               include "allVolVectorEndPoints.H"
                            }
                            else
                            {
#                               include "volVectorEndPoints.H"
                            }
                            break;
                        default:
                            FatalErrorIn("gmtFoam")
                            << fieldsToPlot[ifield].plotType()
                            << " is not a valid plot type for volVectorField "
                            << fieldsToPlot[ifield].name()
                            << ". Valid types are vectors"
                            << exit(FatalError);
                    }
                }
                else if (fieldHeader.headerClassName()=="vectorField")
                {
                    vectorIOField U(fieldHeader);

                    //scalar vectorScale = fieldsToPlot[ifield].vectorScale();
                    label vectorFreq   = fieldsToPlot[ifield].vectorFreq();

                    if
                    (
                        fieldsToPlot[ifield].plotType()
                     == FieldToPlot::VECTOR_END_POINTS
                    )
                    {
#                           include "vectorEndPoints.H"
                    }
                    else
                    {
                            FatalErrorIn("gmtFoam")
                            << fieldsToPlot[ifield].plotType()
                            << " is not a valid plot type for vectorField "
                            << fieldsToPlot[ifield].name()
                            << ". Valid types are vectorEndPoints"
                            << exit(FatalError);
                    }
                }
                else if(fieldHeader.headerClassName()=="surfaceVectorField")
                {
                    surfaceVectorField U(fieldHeader, mesh);

                    if (fieldsToPlot[ifield].plotType() == FieldToPlot::VECTORS)
                    {
                        scalar vectorScale = fieldsToPlot[ifield].vectorScale();
                        label vectorFreq   = fieldsToPlot[ifield].vectorFreq();
                        if (plotAllCells)
                        {
                            #include "allSurfaceVectors.H"
                        }
                        else
                        {
                            #include "surfaceVectors.H"
                        }
                    }
                    else if
                    (
                     fieldsToPlot[ifield].plotType()==FieldToPlot::VECTOR_CONTOURS
                    )
                    {
#                       include "surfaceVectorContours.H"
                    }
                    else if
                    (
                     fieldsToPlot[ifield].plotType()==FieldToPlot::VECTOR_END_POINTS
                    )
                    {
                        scalar vectorScale = fieldsToPlot[ifield].vectorScale();
                        label vectorFreq   = fieldsToPlot[ifield].vectorFreq();
#                       include "surfaceVectorEndPoints.H"
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                        << fieldsToPlot[ifield].plotType()
                        << " is not a valid plot type for surfaceVectorField "
                        << fieldsToPlot[ifield].name()
                        << ". Valid types are vectors vectorEndPoints and vectorContours"
                        << exit(FatalError);
                    }
                }
                else if(fieldHeader.headerClassName()=="surfaceScalarField")
                {
                    if
                    (
                     fieldsToPlot[ifield].plotType() == FieldToPlot::RAW_VALUES
                    )
                    {
#                       include "rawSurfaceValues.H"
                    }
                    else if
                    (
                     fieldsToPlot[ifield].plotType() == FieldToPlot::RAW_FLUXES
                    )
                    {
#                       include "rawFluxValues.H"
                    }
                    else if
                    (
                         fieldsToPlot[ifield].plotType()
                      == FieldToPlot::FILLED_CONTOURS
                    )
                    {
                        #include "filledSurfaceContours.H"
                    }
                    else if
                    (
                     fieldsToPlot[ifield].plotType() ==FieldToPlot::SOLID_CONTOURS
                    )
                    {
                        scalar colourMin = fieldsToPlot[ifield].min();
                        scalar colourMax = fieldsToPlot[ifield].max();
                        scalar colourStep = fieldsToPlot[ifield].delta();
#                    include "surfaceContours.H"
                    }
                    else if
                    (
                     fieldsToPlot[ifield].plotType() ==FieldToPlot::MESH_RANGE
                    )
                    {
#                    include "meshRange.H"
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                        << fieldsToPlot[ifield].plotType()
                        << " is not a valid plot type for surfaceScalarField "
                        << fieldsToPlot[ifield].name()
                        << ". Valid types are rawValues or contours"
                        << exit(FatalError);
                    }
                }
//                else if(fieldsToPlot[ifield].plotType() == FieldToPlot::ADVECTED_CONTOURS)
//                {
//                    const TRiSKData& triskData = TRiSKData::New(mesh);
//                    const fvMesh& dualMesh = triskData.dualMesh();    
//                    contourList cf(mesh, dualMesh, triskData, 0, dimless);
//                    scalar colourMin = fieldsToPlot[ifield].min();
//                    scalar colourMax = fieldsToPlot[ifield].max();
//                    scalar colourStep = fieldsToPlot[ifield].delta();
//                    {
//#                       include "advectedContours.H"
//                    }
//                }
                else if (fieldsToPlot[ifield].plotType() == FieldToPlot::MESH)
                {
                    if (plotAllCells)
                    {
#                       include "allMeshPlot.H"
                    }
                    else
                    {
#                       include "meshPlot.H"
                    }
                }
                else if(fieldsToPlot[ifield].plotType()==FieldToPlot::MESHPOINTS)
                {
                    if (!plotAllCells)
                    {
#                       include "meshPoints.H"
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                        << " must specify a patchName to use for meshPoints"
                        << exit(FatalError);
                    }
                }
                else if(fieldsToPlot[ifield].plotType()==FieldToPlot::MESHCENTRES)
                {
                    if (plotAllCells)
                    {
#                       include "allMeshCentres.H"
                    }
                    else
                    {
#                       include "meshCentres.H"
                    }
                }
                else if (fieldsToPlot[ifield].plotType()==FieldToPlot::NUMBERED)
                {
                    #include "numberedMeshPlot.H"
                }
                else
                {
                    FatalErrorIn("gmtFoam")
                    << fieldsToPlot[ifield].plotType() << " is not a valid plot type, specify one of filledContours, solidContours, dashedContours, vectors, rawValues or mesh"
                    << exit(FatalError);
                }
            }
        }

        // Add boundary
        if (plotBoundaryDots)
        {
#           include "addBoundary.H"
        }

        // Add annotations
#       include "annotate.H"

        // finishing off
        if (lastPlot)
        {
            systemCall = "psbasemap -R" + region + " -J"
                              + projection
                              + " -B" + boundaryMarks
                              + " -O >>"
                              + epsFile;
            systemVerbose(systemCall);
            
            if (!overlay)
            {
                // convert to pdf and remove ps file
                systemCall = "ps2pdf " + epsFile + " " + pdfFile + ".pdf";
                systemVerbose(systemCall);
                systemCall = "pdfcrop " + pdfFile + ".pdf " + pdfFile;
                systemVerbose(systemCall);
    
                systemCall = "rm " + pdfFile + ".pdf " + epsFile;
                systemVerbose(systemCall);
            }
        }
    }

    Info<< "\n end \n";

    return(0);
}


// *************************************************************************
