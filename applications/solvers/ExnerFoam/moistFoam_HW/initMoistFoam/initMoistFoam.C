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
    initMoistFoam_HW

Description
    with inital uniform theta_e and uniform r_t
    Initialise moistFoam, the compressible atmospheric solver with moisture

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "atmosphere.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readInitProps.H"
    HodgeOps H(mesh);
    surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    static scalar piby2 = 0.5*constant::mathematical::pi;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    bool innerConverged = false;
    bool outerConverged = false;
    topBCval = Exner.boundaryField()[topBC][0];
    for
    (
        label iter = 0;
        iter < maxIters && !(innerConverged && outerConverged);
        iter++
    )
    {
        Info << "Outer iteration " << iter << endl;
        
        // Inner iterations with fixed top boundary
        innerConverged = false;
        for(label BCiter = 0; BCiter < BCiters && !innerConverged; BCiter++)
        {
            Info << "Outer " << iter << " inner " << BCiter << endl;
            gradPcoeff = fvc::interpolate
            (
                air.Cp()*theta*atmos.rhoR()
               /(atmos.sumDensity()*air.R()*atmos.volGas())
            );
            
            fvScalarMatrix ExnerEqn
            (
                fvc::div(un)
              - fvm::laplacian(gradPcoeff, Exner)
            );
            innerConverged
                = ExnerEqn.solve((Exner.name())).nIterations() == 0;

            // More inner iterations for moisture variables
            p = air.pFromExner(Exner);

            for(label it = 0; it < 4; it++)
            {
                // Temperature from thetae0
                atmos.TfromThetae(T, thetae0, rt0);

                rho = atmos.rhoFromP(p,T);

                // Set water vapour to be saturated (with under-relaxation)
                water.gas().rho() == (1-underRelax)*water.gas().rho()
                                  + underRelax*min
                (
                    water.pSat(T)/(T*water.gas().R()),
                    rt0*air.rho()
                );
 
                // Update liquid water
                water.liquid().v() = (rt0*air.rho() - water.gas().rho())
                            /water.liquid().rho();
                
                // Update dry air density
                air.rho() += rho - atmos.sumDensity();
            }
            
            // Theta for hydrostatic balance
            theta == T/Exner;
            
            Info << "T goes from " << min(T).value() << " to " <<max(T).value()
                 << "\nwater vapour tho goes from "
                 << min(water.gas().rho()).value()
                 << " to " << max(water.gas().rho()).value() << endl;
        }
        scalar maxGroundExner = max(Exner.boundaryField()[groundBC]);
        outerConverged = (mag(initExnerGround-maxGroundExner)< BCtol);

        // modify the top boundary condition
        Info << "Exner ground value = " << maxGroundExner
             << "  ground value minus one = "
             << maxGroundExner-1
             << " Exner top BC going from  = " << topBCval << " to ";

        topBCval = (1-boundaryRelaxation)*topBCval
                 + boundaryRelaxation*Exner.boundaryField()[topBC][0]
                       /maxGroundExner;
        topBCval = min(max(topBCval, scalar(0)), scalar(1));
        Info << topBCval << endl;
        Exner.boundaryFieldRef()[topBC] == topBCval;
    }

    // Change the top boundary type to be hydrostaticExner
    ExnerBCs[topBC] = "hydrostaticExner";
    Exner = volScalarField
    (
        IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
        Exner,
        ExnerBCs
    );
    Exner.correctBoundaryConditions();
    Exner.write();
    
    // Reference profiles of theta and atmos for calculating differences
    volScalarField thetaRef
    (
        IOobject("thetaRef", runTime.constant(), mesh), theta
    );
    volScalarField airVapourRhoRef
    (
        IOobject("airVapourRhoRef", runTime.constant(), mesh),
        air.rho()
    );
    volScalarField waterVapourRhoRef
    (
        IOobject("waterVapourRhoRef", runTime.constant(), mesh),
        water.gas().rho()
    );
    volScalarField waterLiquidFracRef
    (
        IOobject("waterLiquidFracRef", runTime.constant(), mesh),
        water.liquid().v()
    );
    thetaRef.write();
    airVapourRhoRef.write();
    waterVapourRhoRef.write();
    waterLiquidFracRef.write();
    
    Info << "Adding bouyancy perturbation and re-calculating water variables"
         << endl;
    // Cell where bouancy perturbation is maximum
    label cellMax = 0;
    forAll(theta, celli)
    {
        scalar L = Foam::sqrt
        (
            sqr((mesh.C()[celli].x() - bubbleCentre.x())/bubbleRadii.x())
          + sqr((mesh.C()[celli].y() - bubbleCentre.y())/bubbleRadii.y())
          + sqr((mesh.C()[celli].z() - bubbleCentre.z())/bubbleRadii.z())
        );

        if (L < 1)
        {
            thetaScale[celli] += thetaPrime*sqr(Foam::cos(piby2*L));
            if (thetaScale[celli] > thetaScale[cellMax])
            {
                cellMax = celli;
            }
        }
    }
    volScalarField thetaRho = theta*atmos.rhoR()
                            /(atmos.sumDensity()*air.R()*atmos.volGas());
    thetaRho*= thetaScale;
    // Iterate so that Exner, and water variables are in balance
    for(label it = 0; it < 50; it++)
    {
        theta = thetaRho/atmos.rhoR()*atmos.sumDensity()*air.R()*atmos.volGas();
        T == theta*Exner;
        rho = atmos.rhoFromP(p,T);

        // Set water vapour to be saturated
        water.gas().rho() == min
        (
            water.pSat(T)/(T*water.gas().R()),
            rt0*air.rho()
        );

        // Update liquid water
        water.liquid().v() = (rt0*air.rho() - water.gas().rho())
                    /water.liquid().rho();
        
        air.rho() += rho - atmos.sumDensity();
    }

    theta.write();
    atmos.write(); 
    
    return 0;
}


// ************************************************************************* //
