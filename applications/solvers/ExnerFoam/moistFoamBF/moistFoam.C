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
#include "moistThermo.H"
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
    #include "createFields.H"
    #include "initContinuityErrs.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
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
            #include "moisture.H"
            #include "rhoThetaEqn.H"
            #include "exnerEqn.H"
            p = pFromExner(Exner, kappa, pRef);
            T = TFromTheta(theta, Exner);
        }
        #include "rhoEqn.H"
        #include "rhoThetaEqn.H"

        // Updates for next time step
        dVdt += rhof*gd
             - H.magd()*rhof*gradPCoeff(theta, rv, rl, Cp, epsilon)
                *fvc::snGrad(Exner)
             - muSponge*V;
        
        #include "compressibleContinuityErrs.H"

        thetae = thetaeFromPrimitive(T,p,Lv,rv,rl,pRef,epsilon,Cp,Cpl,R);
        Info << "rl goes from " << min(rl.internalField())
             << " to " << max(rl.internalField()) << endl;
        Info << "rv goes from " << min(rv.internalField())
             << " to " << max(rv.internalField()) << endl;
        Info << "thetae goes from " << min(thetae.internalField())
             << " to "<<max(thetae.internalField()) << endl;
        Info << "condenseRate goes from "
             << min(condenseRate.internalField())
             << " to " << max(condenseRate.internalField()) << endl;

       runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
