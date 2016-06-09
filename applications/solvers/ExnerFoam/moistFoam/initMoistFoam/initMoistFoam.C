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
            Rm = R + Rv*rv;
            cpml = Cp + Cpv*rv;
            kappam = Rm/cpml;
            gradPcoeff = R*fvc::interpolate(thetaRho/kappam);
            fvScalarMatrix ExnerEqn
            (
                fvc::div(U)
              - fvm::laplacian(gradPcoeff, Exner)
            );
            innerConverged
                = ExnerEqn.solve(mesh.solver(Exner.name())).nIterations() == 0;

            // More inner iterations for moisture variables
            p = pRef*pow(Exner, 1/kappam);
            for(label it = 0; it < 4; it++)
            {
                // Temperature from thetae0
                T = thetae0/pow
                    (
                        p/pRef*epsilon/(rv + epsilon),
                       -R/(Cp+Cpl*(rv+rl))
                    )
                    / Foam::exp(Lv*rv/((Cp+Cpl*(rv+rl))*T));

                // Latent heat and saturation vapour pressure from T
                Lv = Lv0 - (Cpl - Cpv)*(T - T0);
                es = Pcc*pRef*Foam::exp(-Lv0/Rv*(1/T - 1/T0));
                volScalarField qvs = epsilon*es/p
                                    *(1 + qv*(1-epsilon)/epsilon - ql);
 
                // Update qv with under-relaxation
                qv == max
                (
                    min((1-underRelax)*qv + underRelax*qvs, qt0),
                    scalar(0)
                );
                rv == qv/(1-qt0);
                rl == rt0 - rv;
                Info << "Inner iteration, T goes from " << min(T).value()
                     << " to " << max(T).value() << nl;
            }
            ql == qt0 - qv;
        
            // Calculate thetaRho for hydrostatic balance
            thetaRho == T/Exner*(1 + rv/epsilon)/(1+rt0);
            
            Info << "T goes from " << min(T).value() << " to " <<max(T).value()
                 << "\nthetaRho goes from " << min(thetaRho).value()
                 << " to " << max(thetaRho).value() << endl;
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
        
        Info << "ql from " << min(ql).value() << " to " << max(ql).value() 
             << " qv from " << min(qv).value() << " to " << max(qv).value()
             << endl;
    }

    // Reference profiles of thetaRho, qv and ql for calculating differences
    volScalarField thetaRhoRef
    (
        IOobject("thetaRhoRef", runTime.constant(), mesh), thetaRho
    );
    volScalarField qvRef(IOobject("qvRef", runTime.constant(), mesh), qv);
    volScalarField qlRef(IOobject("qlRef", runTime.constant(), mesh), ql);
    thetaRhoRef.write();
    qvRef.write();
    qlRef.write();

    Info << "Adding bouyancy perturbation and re-calculating qv and ql" << endl;
    // Cell where bouancy perturbation is maximum
    label cellMax = 0;
    forAll(thetaRho, celli)
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
//    // Interate so that new thetaRho, qv and ql are in balance
//    thetaRho *= thetaScale;
//    for(label iter = 0; iter < 5; iter++)
//    {
//        T = thetaRho*Exner/(1 + qv*(Rv/R-1) - ql);
//        es = Pcc*pRef*Foam::exp(-Lv0/Rv*(1/T - 1/T0));
//        qv == max
//        (
//            min((1-underRelax)*qv + underRelax*es/p, qt0),
//            scalar(0)
//        );
//    }
    rv == qv/(1-qt0);
    rl == rt0 - rv;
    ql == max(qt0 - qv, scalar(0));

    thetaRho.write();
    qv.write();
    ql.write();
    rv.write();
    rl.write();
    
    // Additional diagnositcs for de-bugging
    p.write();
    T.write();
    Lv.write();
    es.write();

    // Change the top boundary type to be fixedFluxBuoyantExnerMoist
    wordList ExnerBCs = Exner.boundaryField().types();
    ExnerBCs[topBC] = "fixedFluxBuoyantExnerMoist";
    volScalarField ExnerNew
    (
        IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
        Exner,
        ExnerBCs
    );
    ExnerNew.correctBoundaryConditions();
    ExnerNew.write();
    
    // Calculate thetae
    volScalarField thetae
    (
        "thetae",
        T*pow(p*epsilon/(rv + epsilon)/pRef, -R/(Cp+Cpl*rt0))
       *Foam::exp(Lv*rv/((Cp+Cpl*rt0)*T))
    );
    Info << "thetae goes from " << min(thetae).value() << " to "
         << max(thetae).value() << endl;
    thetae.write();

    return 0;
}


// ************************************************************************* //
