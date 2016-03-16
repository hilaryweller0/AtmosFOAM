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
    initMoistFoam

Description
    with inital uniform theta_e and uniform r_t
    Initialise moistFoam, the compressible atmospheric solver with moisture

\*---------------------------------------------------------------------------*/

#include "Hops.H"
#include "fvCFD.H"
#include "ExnerTheta.H"
#include "mathematicalConstants.H"

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
    #include "readInitProps.H"
    Hops H(mesh);
    surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    static scalar piby2 = 0.5*constant::mathematical::pi;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    bool innerConverged = false;
    bool outerConverged = false;
    scalar topBCval = Exner.boundaryField()[topBC][0];
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
            fvScalarMatrix ExnerEqn
            (
                fvc::div(U)
              - fvm::laplacian(Cp*thetaRho, Exner)
            );
            innerConverged
                = ExnerEqn.solve(mesh.solver(Exner.name())).nIterations() == 0;

            // More inner iterations for moisture variables
            p = pRef*pow(Exner, 1/kappa);
            for(label it = 0; it < 4; it++)
            {
                // Rosenbrock step to find temperatures, and vapour pressures
                // Partial pressure of dry air
                volScalarField pd = p*epsilon/(rv + epsilon);

                // Temperature from thetae using Rosenbrock step
                volScalarField a = pow(pd/pRef, -R/(Cp+Cpl*rt0));
                volScalarField b = Lv*rv/(Cp + Cpl*rt0);
                T -= (T*a - thetae*Foam::exp(-b/T))/(a*(1 - b/T));
                Lv = Lv0 - (Cpl - Cpv)*(T - T0);
        
                // Calculate saturation vapour pressure and set rv to be saturated
                es = Pcc*pRef*Foam::exp(TccScale*(T - T0)/(T + Tcc));
                // Update rv with under-relaxation
                rv == (1-underRelax)*rv + underRelax*epsilon*(es/(p-es));
                
                volScalarField rvs = epsilon*(es/(p-es));
                condenseRate = (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T))));
                Info << "condenseRate goes from "  << min(condenseRate).value() 
                     << " to " << max(condenseRate).value() << endl;
            }
            rl == max(rt0 - rv, scalar(0));
        
            // Calculate theta and thetaRho for hydrostatic balance
            theta == T/Exner;
            thetaRho = fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl),"theta");
        }
        scalar maxGroundExner = max(Exner.boundaryField()[groundBC]);
        outerConverged = (mag(1-maxGroundExner)< BCtol);

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
        Exner.boundaryField()[topBC] == topBCval;
        
        Info << "thata goes from " << min(theta).value()
             << " to " << max(theta).value() << " rl from "
             << min(rl).value() << " to " << max(rl).value() 
             << " rv from " << min(rv).value() << " to " << max(rv).value()
             << endl;
    }

    // Reference profiles of theta, rv and rl for calculating differences
    volScalarField thetaRef
    (
        IOobject("thetaRef", runTime.constant(), mesh), theta
    );
    volScalarField rvRef(IOobject("rvRef", runTime.constant(), mesh), rv);
    volScalarField rlRef(IOobject("rlRef", runTime.constant(), mesh), rl);
    thetaRef.write();
    rvRef.write();
    rlRef.write();

    Info << "Adding bouyancy perturbation and re-calculating rv and rl" << endl;
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
            theta[celli] *= 1 + (thetaPrime*sqr(Foam::cos(piby2*L))/T0).value();
        }
    }
    T = theta*Exner;
    Lv = Lv0 - (Cpl - Cpv)*(T - T0);
    es = Pcc*pRef*Foam::exp(TccScale*(T - T0)/(T + Tcc));
    rv == epsilon*(es/(p-es));
    rl == max(rt0 - rv, scalar(0));

    theta.write();
    rv.write();
    rl.write();
    
//    // Additional diagnositcs for de-bugging
//    p.write();
//    T.write();
//    Lv.write();
//    es.write();

    // Change the top boundary type to be fixedFluxBuoyantExner
    wordList ExnerBCs = Exner.boundaryField().types();
    ExnerBCs[topBC] = "fixedFluxBuoyantExner";
    volScalarField ExnerNew
    (
        IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
        Exner,
        ExnerBCs
    );
    ExnerNew.correctBoundaryConditions();
    ExnerNew.write();

    return 0;
}


// ************************************************************************* //
