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
    partitionedMoistFoam

Description
    Transient Solver for buoyant, inviscid, compressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and V (flux in d direction)
    and moisture with all variables partitioned and solved using 
    conditionally averaged equations

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "partitionedAtmosphere.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    HodgeOps H(mesh);
    surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initContinuityErrs.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
//    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
//    const int nNonOrthCorr =
//        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const int nThetaCorr = itsDict.lookupOrDefault<int>("nThetaCorr", 2);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "phaseEqns.H"
            for(int thetaCorr = 0; thetaCorr < nThetaCorr; thetaCorr++)
            {
                #include "rhoThetaEqn.H"
            }
            atmosParts.updateSigmas(p);
            //#include "exnerEqn.H"
            p = air.pFromExner(Exner);
        }

//        #include "compressibleContinuityErrs.H"

       runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
