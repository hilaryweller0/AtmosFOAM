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
    moistFoamBF

Description
    Transient Solver for buoyant, inviscid, incompressible, non-hydrostatic flow
    using the advective equations as described by Brian and Fritsch.

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
    #include "readThermoPropertiesMoist.H"
    Hops H(mesh);
    surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    dimensionedScalar acousticCo("acousticCo", speedSound*0.5*dt*max
    (
        fvc::surfaceSum(mesh.magSf())()/mesh.V()
    ));
    Info << "Maximum accoustic Courant number = " << acousticCo << nl << endl;
    
    #include "createFields.H"
    const dimensionedScalar initHeat = fvc::domainIntegrate(theta*rho);
    #include "initEnergy.H"
    #include "energy.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
//    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
//    const int nNonOrthCorr =
//        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        // update old time variables for Crank-Nicholson
        V.oldTime() += (1-offCentre)*dt*dVdt;

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "rhoEqn.H"
            #include "rhoThetaEqn.H"
            #include "moisture.H"
            // Equation of state
            p = rho*R*T*(1 + rv/epsilon);
            Exner = pow(p/pRef, kappa);
            #include "momentum.H"
        }
        
        // Updates for next time step
        
        thetaRho = fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl), "theta");
        divU = fvc::div(U);
        divUtheta = fvc::div(U, theta);
        divUrv = fvc::div(U, rv);
        divUrl = fvc::div(U, rl);

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho)
                                        - initHeat;
        normalisedHeatDiff = (totalHeatDiff/initHeat).value();
        #include "energy.H"

        thetae = T*pow(p/pRef*epsilon/(rv + epsilon), -R/(Cp+Cpl*(rv+rl)))
                *Foam::exp(Lv*rv/((Cp+Cpl*(rv+rl))*T));

        Info << "rl goes from " << min(rl).value() << " to "<<max(rl).value()
             <<endl;
        Info << "rv goes from " << min(rv).value() << " to "<<max(rv).value()
             << endl;

       runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
