/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    exnerFoamH

Description
    Transient Solver for buoyant, inviscid, incompressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and V (flux in d direction)

\*---------------------------------------------------------------------------*/
#include "HodgeOps.H"
#include "fvCFD.H"
#include "ExnerTheta.H"
#include "OFstream.H"
#include "gravity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createGravity.H"
    #include "readRotation.H"
    #include "readThermoProperties.H"
    HodgeOps H(mesh);
    surfaceScalarField gd("gd", g() & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initContinuityErrs.H"
    const dimensionedScalar initHeat = fvc::domainIntegrate(theta*rho);
    #include "initEnergy.H"
    #include "energy.H"
    #include "initCourantFile.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
<<<<<<< HEAD:applications/solvers/ExnerFoam/ExnerFoamCPinterpGrad/exnerFoam.C
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const dimensionedScalar radiativeTimescale("radiativeTimescale", dimTime, 60*60*24);
=======
    const int nThetaCorr = itsDict.lookupOrDefault<int>("nThetaCorr", 2);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
>>>>>>> partitionedFoam:applications/solvers/ExnerFoam/partitionedMoistFoam/partitionedMoistFoam.C

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
<<<<<<< HEAD:applications/solvers/ExnerFoam/ExnerFoamCPinterpGrad/exnerFoam.C
            #include "rhoEqn.H"
            #include "exnerEqn.H"
        }
        
        #include "rhoEqn.H"
        {
            thetaf += dt * (radiationf - thetaf)/radiativeTimescale;
            bf = thetaf * g.unitFaceNormal();
            b = fvc::reconstruct(bf * mesh.magSf());
            theta == (b & g.unit());
            thetaf = mag(bf) + (1.0 - mag(g.unitFaceNormal()))*fvc::interpolate(theta, "thetaFromb");
        }
        #include "compressibleContinuityErrs.H"
=======
            #include "phaseEqns.H"
            for(int thetaCorr = 0; thetaCorr < nThetaCorr; thetaCorr++)
            {
                #include "rhoThetaEqn.H"
            }
            #include "exnerEqn.H"
            p = air.pFromExner(Exner);
            //atmosParts.updateSigmas(p);
        }

//        #include "compressibleContinuityErrs.H"
>>>>>>> partitionedFoam:applications/solvers/ExnerFoam/partitionedMoistFoam/partitionedMoistFoam.C

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho)
                                        - initHeat;
        normalisedHeatDiff = (totalHeatDiff/initHeat).value();
        Info << "Heat error = " << normalisedHeatDiff << endl;
        #include "energy.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
