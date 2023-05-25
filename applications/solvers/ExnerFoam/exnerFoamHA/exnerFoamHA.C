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
    exnerFoamA

Description
    Transient solver for buoyant, viscous, compressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and phi. 
    Separate solutions for components of the velocity.
    Optional turbulence modelling.
    Optional implicit gravity waves and implicit advection.
    Separate momentum equation for w

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
#include "HodgeOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    HodgeOps Hodge(mesh);
    surfaceScalarField gd("gd", g & Hodge.delta());
    #include "readThermo.H"
    #include "createFields.H"
    //#include "initContinuityErrs.H"
    #include "initEnergy.H"
    #include "energy.H"
    
    
    const Switch SIgravityWaves(mesh.schemes().lookup("SIgravityWaves"));
    const Switch implicitU(mesh.schemes().lookup("implicitU"));
    const Switch implicitT(mesh.schemes().lookup("implicitT"));
    const Switch hydrostatic(mesh.schemes().lookup("hydrostatic"));

    const dictionary& itsDict = mesh.solution().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

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

    // Testing
    volVectorField Urhs
    (
        "Urhs",
        rho*fvc::reconstruct(gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf())
        //g - Cp*theta*fvc::grad(Exner)
    );
    Urhs.write();
    volVectorField gradExner("gradExner", fvc::grad(Exner));
    gradExner.write();
    
    dimensionedScalar H("H", dimLength, 10e3);
    volScalarField decay
    (
        "decay",
        //exp(-mesh.C().component(2)/H)
        sqr((H - mesh.C().component(2))/H)
    );
    decay.write();
    
    surfaceScalarField HsnGradA
    (
        "HsnGradA",
        (2*mesh.Cf().component(2)/H - 2)*(kdir & mesh.Sf())/mesh.magSf()
    );
    surfaceScalarField HsnGrad("HsnGrad", H*fvc::snGrad(decay) - HsnGradA);
    HsnGrad.write();
    
    volVectorField gradDecayA("gradDecayA", 0*mesh.C()/H);
    gradDecayA.replace(2, 2*mesh.C().component(2)/H - 2);
    gradDecayA.write();
    
    volVectorField gradDecay1("gradDecay1", H*fvc::grad(decay)-gradDecayA);
    gradDecay1.write();
    
    volVectorField gradDecay2
    (
        "gradDecay2",
        H*fvc::reconstruct(fvc::snGrad(decay)*mesh.magSf())-gradDecayA
    );
    gradDecay2.write();
    
    volVectorField gradDecay3
    (
        "gradDecay3",
        H*Hodge.reconstructd
        (
            fvc::snGrad(decay, "uncorr")*Hodge.magd()
        )
        -gradDecayA
    );
    gradDecay3.write();

    //FatalErrorIn("exnerFoamA") << exit(FatalError);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "rhoThetaEqn.H"
            #include "UEqn.H"
            // Exner and momentum equations
            #include "exnerEqn.H"
        }
        #include "rhoThetaEqn.H"
        //#include "compressibleContinuityErrs.H"

        #include "thermoUpdate.H"
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
