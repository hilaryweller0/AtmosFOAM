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
    partitionedExnerFoamAdv_newCN

Description
    Transient Solver for dry, buoyant, compressible, non-hydrostatic
    partitioned flow, advective form momentum equation, with optional turbulence
    modelling.
    Uses OpenFOAM's in-built functionality for Crank-Nicholson timestepping 
    rather than the top-level timestepping usually utilised in AtmosFOAM, hence
    "newCN".

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "ExnerTheta.H"
#include "rhoThermo.H"
#include "Partitioned.H"
#include "PartitionedFields.H"
#include "CrankNicolsonDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    #include "readTransferCoeffs.H"
    HodgeOps H(mesh);
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initDiags.H"
    #include "calcDiags.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    fv::CrankNicolsonDdtScheme<vector> drhoUdt
    (
        mesh,
        mesh.schemesDict().subDict("ddtSchemes").lookup("volFlux_CN")
    );
    const scalar ocCoeff = drhoUdt.ocCoeff();

    // Validate turbulence fields after construction and update derived fields
    for(label ip = 0; ip < nParts; ip++)
    {
        turbulence[ip].validate();
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr+1; ucorr++)
        {
            #include "rhoSigmaEqn.H"
            //#include "massTransfers.H"
            #include "thetaEqn.H"
            #include "sigma.H"
            //#include "calculateDrag.H"
            if (ucorr < nOuterCorr)
            {
                #include "exnerEqn.H"
            }
        }
        
        //- Solve the turbulence equations and correct the turbulence viscosity
        for(label ip = 0; ip < nParts; ip++)
        {
            turbulence[ip].correct();
        }

        #include "compressibleContinuityErrs.H"
        //#include "correctContinuityErrs.H"
        #include "calcDiags.H"
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
