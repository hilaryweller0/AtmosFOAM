{
    volScalarField rhoSum = rho[0];
    for(label ip = 1; ip < nParts; ip++)
    {
        rhoSum += rho[ip];
    }
    
    thetaSum = sigma[0]*pow
    (
        rho[0]*R*theta[0]/pRef, kappa/(1-kappa)
    );
    for(label ip = 1; ip < nParts; ip++)
    {
        thetaSum += sigma[ip]*pow
        (
            rho[ip]*R*theta[ip]/pRef, kappa/(1-kappa)
        );
    }
    thetaMean = pRef/(rhoSum*R)*pow(thetaSum, (1-kappa)/kappa);

    Exner = pow(rhoSum*R*thetaMean/pRef, kappa/(1-kappa));
    PsiSum = rhoSum/Exner;

    gradPcoeff[0] = Cp*fvc::interpolate(sigmaRho[0])*fvc::interpolate(theta[0]);
    gradPcoeffSum = gradPcoeff[0];
    for(label ip = 1; ip < nParts; ip++)
    {
        gradPcoeff[ip] += Cp*fvc::interpolate(sigmaRho[ip])
                            *fvc::interpolate(theta[ip]);
        gradPcoeffSum += gradPcoeff[ip];
        Psi[ip] = rho[ip]/Exner;
    }

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition
        fluxSum
             == dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0));
        for(label ip = 0; ip < nParts; ip++)
        {
            // Rate of change of flux without gravity
            dFluxdt[ip] = -mesh.Sf() &
            (
                fvc::interpolate(fvc::div(flux[ip], u[ip]))
              + (2*fvc::interpolate(sigmaRho[ip])*(Omega ^ Uf[ip]))
            );
            // Flux with gravity
            flux[ip] = flux[ip].oldTime()
                     + dt*dFluxdt[ip]
                     - dt*fvc::interpolate(sigmaRho[ip])*gSf;

            fluxSum += flux[ip];
        }
        // Needed for Exner boundary conditions
        surfaceScalarField rhog("rhog", fvc::interpolate(rhoSum)*gSf);

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(PsiSum, Exner)
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
            if (nonOrth == nNonOrthCorr)
            {for(label ip = 0; ip < nParts; ip++)
            {
                flux[ip] -= dt*gradPcoeff[ip]*fvc::snGrad(Exner)*mesh.magSf();
                rho[ip] = Psi[ip]*Exner;
                u[ip] = fvc::reconstruct
                (
                    flux[ip]/fvc::interpolate(sigmaRho[ip])
                );
                Uf[ip] = fvc::interpolate(u[ip]);
                dFluxdt[ip] -= fvc::interpolate(sigmaRho[ip])*gSf
                            - gradPcoeff[ip]*fvc::snGrad(Exner)*mesh.magSf();
            }}
        }
    }
}
