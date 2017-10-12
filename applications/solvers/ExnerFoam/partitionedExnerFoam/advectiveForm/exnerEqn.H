{
    Psi = pow(rhoSum,(2*kappa-1)/(kappa-1))
         *pow(R/pRef*rhoTheta/rhoSum, kappa/(kappa-1));
    
    gradPcoeffSum = Cp*fvc::interpolate(rhoTheta);

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition
        fluxSum
             == dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0));
        for(label ip = 0; ip < nParts; ip++)
        {
            // Rate of change of flux without gravity
            dVolFluxdt[ip] = (mesh.Sf() &
            (
              -fvc::interpolate
               (
                   fvc::div(volFlux[ip], u[ip]) - u[ip] * fvc::div(volFlux[ip]),
                   "convection"
               )
             - 2*(Omega ^ Uf[ip])
            ));
            // Flux with gravity
            volFlux[ip] = volFlux[ip].oldTime() + dt*dVolFluxdt[ip] + dt*gSf;
            flux[ip] = fvc::interpolate(sigmaRho[ip])*volFlux[ip];
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

            if(ucorr==nOuterCorr-1 && nonOrth==nNonOrthCorr && corr==nCorr-1)
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
                volFlux[ip] -= dt*Cp*fvc::interpolate(theta[ip])
                               *fvc::snGrad(Exner)*mesh.magSf();
                flux[ip] = fvc::interpolate(sigmaRho[ip])*volFlux[ip];
                u[ip] = fvc::reconstruct(volFlux[ip]);
                Uf[ip] = fvc::interpolate(u[ip]);
                dVolFluxdt[ip] += gSf - Cp*fvc::interpolate(theta[ip])
                                        *fvc::snGrad(Exner)*mesh.magSf();
            }}
        }
    }
}
