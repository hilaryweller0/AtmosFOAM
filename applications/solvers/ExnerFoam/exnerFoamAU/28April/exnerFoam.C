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
    Need to include implicit gravity waves and implicit advection

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "physicalProperties.H"
#include "fundamentalConstants.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "constTransport.H"
#include "OFstream.H"
#include "rhoThermo.H"
#include "EulerDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    //#include "createSponge.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermo.H"
    #include "createFields.H"
    //#include "initContinuityErrs.H"
    const dimensionedScalar initHeat = fvc::domainIntegrate(theta*rho);
    #include "initEnergy.H"
    #include "energy.H"
    
    const Switch SIgravityWaves(mesh.schemes().lookup("SIgravityWaves"));
    const Switch impU(mesh.schemes().lookup("implicitU"));
    const Switch stagger(mesh.schemes().lookup("stagger"));
    const dictionary& itsDict = mesh.solution().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const Switch hydrostatic(mesh.schemes().lookup("hydrostatic"));
    const scalar ocCoeff
    (
        readScalar(mesh.schemes().subDict("ddtSchemes").lookup("ocCoeff"))
    );
    const scalar ocAlpha = 1/(1+ocCoeff);
    // Pre-defined time stepping scheme
    fv::EulerDdtScheme<scalar> EulerDdt(mesh);
    
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
            #include "UEqn.H"
            // Exner and momentum equations
            #include "exnerEqn.H"
        }
        #include "rhoThetaEqn.H"
        //#include "compressibleContinuityErrs.H"

        // Update the pressure and temperature based on the new Exner
        thermo.p() = pRef*pow(Exner, 1/kappa);
        thermo.T() = theta*Exner;
        thermo.he() == thermo.he(thermo.p(),thermo.T());
        thermo.correct();

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho) - initHeat;
        Info << "Heat error = " << (totalHeatDiff/initHeat).value() << endl;
        #include "energy.H"
        
        //- Solve the turbulence equations and correct the turbulence viscosity
        turbulence->correct(); 

        if (runTime.writeTime())
        {
            runTime.write();
            surfaceVectorField Uf("Uf", linearInterpolate(U));
            Uf += (phi/rhof - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
            Uf.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
