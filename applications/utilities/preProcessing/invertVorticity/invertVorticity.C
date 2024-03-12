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
    invertVorticity.C

Description
    Reads in the scalar valued vorticity on a 2d mesh and inverts to find the
    streamfunciton and velocity. Needs to read in the velocity field, Uf for
    the boundary conditions

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvScalarMatrix.H"
#include "fvmLaplacian.H"
#include "fvcCurl.H"
#include "fvcCurlf.H"
#include "fvcVolumeIntegrate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("dictionary name (in system)");
    timeSelector::addOptions();
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
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
            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word dictName = args.args()[1].c_str();
    Info << "Read in constant background wind from" << dictName << endl;
    IOdictionary initDict
    (
        IOobject
        (
            dictName, mesh.time().system(), mesh, IOobject::MUST_READ
        )
    );
    const dimensionedVector U0(initDict.lookup("U0"));
    const Switch setReference(initDict.lookup("setStreamFunctionReference"));

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        Info << "Mesh has normal direction" << flush;
        const vector meshNormal = 0.5*(Vector<label>(1,1,1)-mesh.geometricD());
        Info << meshNormal << endl;

        Info << "Reading in vorticity" << endl;
        volScalarField vorticity
        (
            IOobject("vorticity", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );

        Info << "Reading the streamFunction" << endl;
        volScalarField streamFunction
        (
            IOobject("streamFunction", runTime.timeName(), mesh,
                     IOobject::MUST_READ),
            mesh
        );
                
        // Set the streamfunction for the background uniform flow
        const dimensionedVector velocityPerp = meshNormal ^ U0;
        if (magSqr(velocityPerp).value() > SMALL)
        {
            streamFunction == magSqr(U0)/magSqr(velocityPerp)
                             *(mesh.C() & velocityPerp);
        }

        // Invert the streamfunction to find the vorticity
        bool converged = false;
        for(label it = 0; it < 10 && !converged; it++)
        {
            fvScalarMatrix streamFuncEqn
            (
                fvm::laplacian(streamFunction) == -vorticity
            );
            if (setReference)
            {
                streamFuncEqn.setReference(0,0);
                // Ensure that the vorticity sums to zero
                dimensionedScalar vorticityMean
                    = fvc::domainIntegrate(vorticity)/gSum(mesh.V());
                Info << "Subtracting vorticityMean = " << vorticityMean << endl;
                streamFuncEqn -= vorticityMean;
            }
            solverPerformance sp = streamFuncEqn.solve();
            converged = sp.nIterations() <= 0;
        }

        streamFunction.write();

        // Set the velocity from the streamfunction
        volVectorField U
        (
            IOobject("U", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
            fvc::curl(streamFunction*meshNormal)
        );
        U == fvc::curl(streamFunction*meshNormal);
        surfaceVectorField Uf
        (
            IOobject("Uf", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
            linearInterpolate(U)
        );
        Uf = linearInterpolate(U);
        
        U.write();
        Uf.write();
        
//        Info << "Projecting the velocity field to be divergence free" << endl;
//        surfaceScalarField phiv("phiv", Uf & mesh.Sf());
//        
//        Info << "Reading the velocity potential" << endl;
//        volScalarField velPot
//        (
//            IOobject("velPot", runTime.timeName(), mesh,
//                     IOobject::MUST_READ),
//            mesh
//        );

//        converged = false;
//        for(label it = 0; it < 10 && !converged; it++)
//        {
//            fvScalarMatrix velPotEqn
//            (
//                fvc::div(phiv) - fvm::laplacian(velPot)
//            );
//            if (setReference)
//            {
//                velPotEqn.setReference(0,0);
//            }
//            solverPerformance sp = velPotEqn.solve();
//            converged = sp.nIterations() <= 0;
//            phiv += velPotEqn.flux();
//        }
//        
//        Uf += (phiv - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
//        Uf.write();
//        velPot.write();
//        phiv.write();
//        volScalarField divU("divU", fvc::div(phiv));
//        divU.write();
    }
    
    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
