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

#include "Hops.H"
#include "fvCFD.H"
#include "ExnerTheta.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "orthogonalBoundaries.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    Hops H(mesh);
    surfaceScalarField gd("gd", g & H.delta());
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
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    const Switch SIgravityWaves(mesh.schemesDict().lookup("SIgravityWaves"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        // update old time variables for Crank-Nicholson
        V.oldTime() += (1-offCentre)*dt*dVdt;
        // Old part of theta change (before any variables are updated)
        if (SIgravityWaves)
        {
            thetaf.oldTime() = fvc::interpolate
            (
                (
                    rho.oldTime()*theta.oldTime()
                  - (1-offCentre)*dt*divUtheta.oldTime()
                )/(rho.oldTime() - (1-offCentre)*dt*divU),
                "interpolate(theta)"
            );
        }

        // u predictor
        //u = H.reconstructd((V.oldTime() + offCentre*dt*dVdt)/rhof);

        for (int ucorr=0; ucorr<nOuterCorr; ucorr++)
        {
            #include "rhoThetaEqn.H"

            // Exner and momentum equations
            #include "exnerEqn.H"
        }
        
        #include "rhoThetaEqn.H"
        
        // Updates for next time step
//        if (Charney-Phillips)
//        {
//            thetaf -= dt*offCentre*V/(rhof*H.magd())*fvc::snGrad(theta);
//            theta = fvc::faceToCell(thetaf);
//        }
//        else
        {
            thetaf = fvc::interpolate(theta);
        }
        dVdt += rhof*gd - H.magd()*Cp*rhof*thetaf*fvc::snGrad(Exner)
              - muSponge*V;
        divU = fvc::div(U);
        divUtheta = fvc::div(U, theta);
        
        #include "compressibleContinuityErrs.H"

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho)
                                        - initHeat;
        Info << "Heat error = " << (totalHeatDiff/initHeat).value() << endl;
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
