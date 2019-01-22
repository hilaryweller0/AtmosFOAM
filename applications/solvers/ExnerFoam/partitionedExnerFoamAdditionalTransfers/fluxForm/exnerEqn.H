{
    Psi = pow(rhoSum,(2*kappa-1)/(kappa-1))
         *pow(R/pRef*rhoTheta/rhoSum, kappa/(kappa-1));
    
//    gradPcoeff[0] = Cp*fvc::interpolate(sigmaRho[0])*fvc::interpolate(theta[0]);
//    gradPcoeffSum = gradPcoeff[0];
//    for(label ip = 1; ip < nParts; ip++)
//    {
//        gradPcoeff[ip] = Cp*fvc::interpolate(sigmaRho[ip])*fvc::interpolate(theta[ip]);
//        gradPcoeffSum += gradPcoeff[ip];
//    }

    surfaceScalarField gradExner
    (
        IOobject("gradExner", runTime.timeName(), mesh), -fvc::snGrad(Exner)
    );

    gradPcoeffSum = Cp*fvc::interpolate(rhoTheta);
//    for (label ip = 0; ip < nParts; ip++)
//    {
//        gradPcoeff[ip] = fvc::interpolate(sigma[ip])*gradPcoeffSum;
//    }

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition
        fluxSum
             == dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0));
        for(label ip = 0; ip < nParts; ip++)
        {
            // Rate of change of flux without gravity
            dFluxdt[ip] = drag[ip] + (mesh.Sf() &
            (
              - fvc::interpolate(fvc::div(flux[ip], u[ip]))
              - (2*fvc::interpolate(sigmaRho[ip])*(Omega ^ Uf[ip]))
            ));
            // Flux with gravity
            flux[ip] = flux[ip].oldTime()
                     + dt*dFluxdt[ip]
                     + dt*fvc::interpolate(sigmaRho[ip])*gSf;

            // Momentum tranfer with mass tranfer
            for(label jp=0; jp < nParts; jp++) {if (ip != jp)
            {
                flux[ip] += dt*
                (
                    flux[jp]*fvc::interpolate(massTransfer[jp][ip]/sigmaRho[jp])
                  - flux[ip]*fvc::interpolate(massTransfer[ip][jp]/sigmaRho[ip])
                );
            }}
            
            fluxSum += flux[ip];
        }
        // Needed for Exner boundary conditions
        surfaceScalarField rhog("rhog", fvc::interpolate(rhoSum)*gSf);

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            // Explicit version
            //Exner = pow(R*rhoTheta/pRef, kappa/(1-kappa));
            // Semi-implicit version
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + fvc::div(fluxSum)
              - fvm::laplacian(dt*gradPcoeffSum, Exner)
            );

            if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
            {
                ExnerEqn.solve(mesh.solver(Exner.name() + "Final"));
            }
            else
            {
                ExnerEqn.solve(mesh.solver(Exner.name()));
            }

            // Back-substitutions
            if (nonOrth == nNonOrthCorr) { for(label ip = 0; ip < nParts; ip++)
            {
                surfaceScalarField sigmaf = fvc::interpolate(sigma[ip]);
                flux[ip] -= dt*sigmaf*gradPcoeffSum*fvc::snGrad(Exner)*mesh.magSf();
                u[ip] = fvc::reconstruct(flux[ip])/sigmaRho[ip];
                Uf[ip] = fvc::interpolate(u[ip]);
                dFluxdt[ip] += fvc::interpolate(sigmaRho[ip])*gSf
                            - sigmaf*gradPcoeffSum*fvc::snGrad(Exner)*mesh.magSf();
            }}
        }
    }
}
