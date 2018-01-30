{
    Psi = pow(rho.sum(),(2*kappa-1)/(kappa-1))
         *pow(R/pRef*theta.sum()/rho.sum(), kappa/(kappa-1));
    
    gradPcoeffSum = Cp*fvc::interpolate(theta.sum(), "theta");

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition (without mass tranfers)
        for(label ip = 0; ip < nParts; ip++)
        {
            // Rate of change of flux without gravity
            volFlux.ddt()[ip] = drag[ip] + (mesh.Sf() &
            (
              -fvc::interpolate
               (
                   fvc::div(volFlux[ip], u[ip], "u")
                 - u[ip] * fvc::div(volFlux[ip]),
                   "convection"
               )
             - 2*(Omega ^ Uf[ip])
            ));
            
            // Flux with gravity
            volFlux[ip] = volFlux[ip].oldTime() + dt*
            (
                (1-offCentre)*volFlux.ddt()[ip].oldTime()
              + offCentre*(volFlux.ddt()[ip] + gSf)
            );
            flux[ip] = sigmaRhof[ip]*volFlux[ip];
        }
        flux.updateSum();
        // Needed for Exner boundary conditions
        surfaceScalarField rhog
        (
            "rhog",
            fvc::interpolate(rho.sum(), "rho")*gSf
        );

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + (1-offCentre)*fvc::div(flux.sum().oldTime())
              + fvc::div(flux.sum())
              - fvm::laplacian(offCentre*dt*gradPcoeffSum, Exner)
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
            if (nonOrth == nNonOrthCorr)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    surfaceScalarField thetaf = fvc::interpolate
                    (
                        theta[ip], "theta"
                    );
                    volFlux[ip] -= offCentre*dt*Cp*thetaf
                                   *fvc::snGrad(Exner)*mesh.magSf();

                    volFlux.ddt()[ip] += gSf - Cp*thetaf
                                        *fvc::snGrad(Exner)*mesh.magSf();
                }
            
                // Include the mass transfer terms implicitly (for nParts = 2)
                volFlux[0] =
                (
                    (1+transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
                )/(1+transferf.sum());
                volFlux[1] = (volFlux[1] + transferf[0]*volFlux[0])
                            /(1+transferf[0]);

                // Update mass flux and velocities based on the flux
                for(label ip = 0; ip < nParts; ip++)
                {
                    flux[ip] = sigmaRhof[ip]*volFlux[ip];
                    u[ip] = fvc::reconstruct(volFlux[ip]);
                    Uf[ip] = fvc::interpolate(u[ip], "Uf");
                }
                flux.updateSum();
                volFlux.updateSum();
            }
        }
    }
}

