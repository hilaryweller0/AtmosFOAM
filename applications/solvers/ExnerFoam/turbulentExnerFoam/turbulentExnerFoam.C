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
    turbulentExnerFoam

Description
    Transient solver for buoyant, viscous, compressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and phi. 
    Optional turbulence modelling.

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "ExnerTheta.H"
#include "OFstream.H"
#include "rhoThermo.H"
#include "CrankNicolsonDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    HodgeOps H(mesh);
//    const surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initContinuityErrs.H"
    const dimensionedScalar initHeat = fvc::domainIntegrate(theta*rho);
    #include "initEnergy.H"
    #include "energy.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    fv::CrankNicolsonDdtScheme<vector> drhoUdt
    (
        mesh,
        mesh.schemesDict().subDict("ddtSchemes").lookup("ddt(rho,U)_CN")
    );
    const scalar ocCoeff = drhoUdt.ocCoeff();
    
    turbulence->validate();   //- Validate turbulence fields after construction
                            //  and update derived fields as required

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        for (int ucorr=0; ucorr<nOuterCorr; ucorr++)
        {
            #include "rhoThetaEqn.H"
//            #include "UEqn.H"

            // Exner and momentum equations
            #include "exnerEqn.H"
        }
        
        #include "rhoThetaEqn.H"
        
        // Update rates of change for next time step
        dPhidt += rhof*(gSf - mesh.magSf()*Cp*thetaf*fvc::snGrad(Exner));
        
        #include "compressibleContinuityErrs.H"

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho)
                                        - initHeat;
        Info << "Heat error = " << (totalHeatDiff/initHeat).value() << endl;
        #include "energy.H"
        
        //- Solve the turbulence equations and correct the turbulence viscosity
        turbulence->correct(); 

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
