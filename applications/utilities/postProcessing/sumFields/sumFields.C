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
    sumFields

Description
    sum the fields in the argument list using the scaling factors given in the argument list

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("outTime");
    Foam::argList::validArgs.append("outFileName");
    Foam::argList::validOptions.insert("scale0","scale0");
    Foam::argList::validArgs.append("time");
    Foam::argList::validArgs.append("fieldName");
    Foam::argList::validOptions.insert("scale1","scale1");
    Foam::argList::validArgs.append("time");
    Foam::argList::validArgs.append("fieldName");
    Foam::argList::validOptions.insert("pow0","pow0");
    Foam::argList::validOptions.insert("pow1","pow1");
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"

    label nFiles = 2;
    List<fileName> fileNames(nFiles);
    scalarList scales(nFiles, scalar(1));
    scalarList powers(nFiles, scalar(1));
    List<fileName> fileTime(nFiles);

    if (args.options().found("scale0"))
    {
        scales[0] = readScalar(IStringStream(args.options()["scale0"])());
    }

    if (args.options().found("scale1"))
    {
        scales[1] = readScalar(IStringStream(args.options()["scale1"])());
    }

    if (args.options().found("pow0"))
    {
        powers[0] = readScalar(IStringStream(args.options()["pow0"])());
    }

    if (args.options().found("pow1"))
    {
        powers[1] = readScalar(IStringStream(args.options()["pow1"])());
    }

    for (label i = 0; i < nFiles; i++)
    {
        fileTime[i] = args.args()[3+2*i];
        fileNames[i] = args.args()[4+2*i];
    }
    scalar timeOut = readScalar(IStringStream(args.args()[1])());
    fileName fileNameOut = args.args()[2];

#   include "createTime.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

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

    IOdictionary envDict
    (
        IOobject
        (
            "environmentalProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    IOobject fieldHeader(fileNames[0], fileTime[0], mesh, IOobject::MUST_READ);

    if ( !fieldHeader.headerOk())
    {
        FatalErrorIn("sumFields") << "Cannot read " << fileNames[0]
            << " from time directory " << fileTime[0] << exit(FatalError);
    }
    else if (fieldHeader.headerClassName() == "volScalarField")
    {
        volScalarField vfIn(fieldHeader, mesh);

        volScalarField vfOut
        (
            IOobject(fileNameOut, runTime.timeName(), mesh),
            scales[0]*pow(vfIn, powers[0])//,
            //zeroGradientFvPatchScalarField::typeName
        );

        if (mag(scales[1]) > SMALL)
        {
            vfIn = volScalarField
            (
                IOobject(fileNames[1],fileTime[1],mesh,IOobject::MUST_READ),
                mesh
            );

            vfOut += scales[1] * pow(vfIn, powers[1]);
        }

        runTime.setTime(timeOut, 0);
        Info << "Writing " << timeOut << "/" << fileNameOut << endl;
        vfOut.write();

        label imin = 0;
        label imax = 0;
        scalar minVal = vfOut[0];
        scalar maxVal = vfOut[0];
        for(label i = 1; i < vfOut.size(); i++)
        {
            if (vfOut[i] < minVal)
            {
                minVal = vfOut[i];
                imin = i;
            }
            else if (vfOut[i] > maxVal)
            {
                maxVal = vfOut[i];
                imax = i;
            }
        }
        Info << "Min = " << minVal << " in cell " << imin
            << " max = " << maxVal << " in cell " << imax << endl;
    }
    else if (fieldHeader.headerClassName() == "volVectorField")
    {
        if (args.options().found("pow0") || args.options().found("pow1"))
        {
            FatalErrorIn("sumFields") << "cannot take a power of a vector field."
                << "\nvector fields " << fileNames[0] << " and " << fileNames[1]
                << " powers " << powers[0] << " and " << powers[1]
                << exit(FatalError);
        }

        volVectorField vfIn(fieldHeader, mesh);

        volVectorField vfOut = scales[0]*vfIn;
        vfOut.rename(fileNameOut);

        if (mag(scales[1]) > SMALL)
        {
            vfIn = volVectorField
            (
                IOobject(fileNames[1],fileTime[1],mesh,IOobject::MUST_READ),
                mesh
            );

            vfOut += scales[1]*vfIn;
        }

        runTime.setTime(timeOut, 0);
        Info << "Writing " << timeOut << "/" << fileNameOut << endl;
        vfOut.write();

        Info << "Max = " << max(mag(vfOut.internalField()))
             << " min = " << min(mag(vfOut.internalField()))
             << endl;
    }
    else if (fieldHeader.headerClassName() == "surfaceScalarField")
    {
        surfaceScalarField sfIn(fieldHeader, mesh);

        surfaceScalarField sfOut = scales[0]*pow(sfIn, powers[0]);
        sfOut.rename(fileNameOut);

        if (mag(scales[1]) > SMALL)
        {
            sfIn = surfaceScalarField
            (
                IOobject(fileNames[1],fileTime[1],mesh,IOobject::MUST_READ),
                mesh
            );

            sfOut += scales[1] * pow(sfIn, powers[1]);
        }

        runTime.setTime(timeOut, 0);
        Info << "Writing " << timeOut << "/" << fileNameOut << endl;
        sfOut.write();

        Info << "Max = " << max(sfOut.internalField())
             << " Min = " << min(sfOut.internalField()) << endl;
    }
    else if (fieldHeader.headerClassName() == "surfaceVectorField")
    {
        if (args.options().found("pow0") || args.options().found("pow1"))
        {
            FatalErrorIn("sumFields") << "cannot take a power of a vector field."
                << "\nvector fields " << fileNames[0] << " and " << fileNames[1]
                << " powers " << powers[0] << " and " << powers[1]
                << exit(FatalError);
        }

        surfaceVectorField vfIn(fieldHeader, mesh);

        surfaceVectorField vfOut = scales[0]*vfIn;
        vfOut.rename(fileNameOut);

        if (mag(scales[1]) > SMALL)
        {
            vfIn = surfaceVectorField
            (
                IOobject(fileNames[1],fileTime[1],mesh,IOobject::MUST_READ),
                mesh
            );

            vfOut += scales[1] * vfIn;
        }

        runTime.setTime(timeOut, 0);
        Info << "Writing " << timeOut << "/" << fileNameOut << endl;
        vfOut.write();

        Info << "Max = " << max(mag(vfOut.internalField()))
             << " Min = " << min(mag(vfOut.internalField()))
             << endl;
    }
    else
    {
        FatalErrorIn("sumFields") << fileNames[0] << " is of type "
            << fieldHeader.headerClassName()
            << " which is not handled by sumFields.\n"
            << "sumFields only handles volScalarFields, surfaceScalarFields, surfaceVectorFields and volVectorFields"
                    << exit(FatalError);
    }

    Info << endl;

    return(0);
}


// ************************************************************************* //
