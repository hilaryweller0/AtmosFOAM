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
    setScalarOverOrography

Description
    Set the initial scalar, T for scalar transport over orography for any time
    or for all of the times in the case directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "VelocityProfile.H"
#include "Mountain.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "tracerFieldFileName", "filename", 
        "specify the name of the tracer field file name (default 'T')"
    );
#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList Times = runTime.times();
    if (Times.size() == 1 && !args.optionFound("constant"))
    {
        Times.append(instant(scalar(0)));
    }

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "Reading initial conditions" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "setScalarOverOrographyDict", runTime.constant(), mesh,
            IOobject::MUST_READ
        )
    );
    
    IOdictionary velocityDict
    (
        IOobject("velocityFieldDict", runTime.constant(), mesh, IOobject::MUST_READ)
    );
    
    // Initial maximum tracer value
    const scalar rho0(readScalar(initDict.lookup("rho0")));
    // Initial tracer position
    const scalar x0(readScalar(initDict.lookup("x0")));
    const scalar z0(readScalar(initDict.lookup("z0")));
    // Half widths
    const scalar Ax(readScalar(initDict.lookup("Ax")));
    const scalar Az(readScalar(initDict.lookup("Az")));

    const string tracerFieldFileName = args.options().found("tracerFieldFileName") ?
                                       args.options()["tracerFieldFileName"] : "T";

    // Comment from Hilary: why are you using a pointer here?
    VelocityProfile* velocityProfile = VelocityProfile::lookup(velocityDict);

    Info << "Creating initial tracer field " << tracerFieldFileName << endl;
    volScalarField T
    (
        IOobject(tracerFieldFileName, runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(tracerFieldFileName, dimless, scalar(0)),
        "zeroGradient"
    );

    Info << "Creating initial tracer field " << (tracerFieldFileName + "f") << endl;
    surfaceScalarField Tf
    (
        IOobject(tracerFieldFileName + "f", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(tracerFieldFileName + "f", dimless, scalar(0)),
        "fixedValue"
    );
    
    // Set the tracer for each time specified
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;
            // Centre of the tracer for this time step
        const point& advectedPt = velocityProfile->pointAtTime(point(x0, 0, z0), runTime.value());

        // Calculating T
        forAll(T, cellI)
        {
            const point& c = mesh.C()[cellI];
            
            // Define r as used in the initial tracer field
            scalar r = Foam::sqrt(sqr((c.x()-advectedPt.x())/Ax)+sqr((c.z()-advectedPt.z())/Az));
        
            if (r <= 1)
            {
                T[cellI] = rho0*sqr(Foam::cos(M_PI*r/2));
            }
            else T[cellI] = 0;
        }
        T.correctBoundaryConditions();
        T.write();

        // Calculating Tf
        forAll(Tf, faceI)
        {
            const point& c = mesh.Cf()[faceI];
            
            // Define r as used in the initial tracer field
            scalar r = Foam::sqrt(sqr((c.x()-advectedPt.x())/Ax)+sqr((c.z()-advectedPt.z())/Az));
        
            if (r <= 1)
            {
                Tf[faceI] = rho0*sqr(Foam::cos(M_PI*r/2));
            }
            else Tf[faceI] = 0;
        }
        Tf.write();
    }
    
    Info<< "End\n" << endl;
    
    return(0);
}

// ************************************************************************* //
