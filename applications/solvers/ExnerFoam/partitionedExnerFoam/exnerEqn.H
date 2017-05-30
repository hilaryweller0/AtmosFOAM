{
    volScalarField rhoSum = rho[0];
    for(label ip = 1; ip < nParts; ip++)
    {
        rhoSum += rho[ip];
    }
    
    rhoTheta = sigma[0]*rho[0]*theta[0];
    for(label ip = 1; ip < nParts; ip++)
    {
        rhoTheta += sigma[ip]*rho[ip]*theta[ip];
    }

    Exner = pow(R*rhoTheta/pRef, kappa/(1-kappa));
    Psi = rhoSum/Exner;

    gradPcoeff[0] = Cp*fvc::interpolate(sigmaRho[0])
                    *fvc::interpolate(theta[0]);
    gradPcoeffSum = gradPcoeff[0];
    for(label ip = 1; ip < nParts; ip++)
    {
        gradPcoeff[ip] += Cp*fvc::interpolate(sigmaRho[ip])
                            *fvc::interpolate(theta[ip]);
        gradPcoeffSum += gradPcoeff[ip];
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
                     + dt*fvc::interpolate(sigmaRho[ip])*gSf;

            fluxSum += flux[ip];
        }
        // Needed for Exner boundary conditions
        surfaceScalarField rhog("rhog", fvc::interpolate(rhoSum)*gSf);

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
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
                flux[ip] -= dt*gradPcoeff[ip]*fvc::snGrad(Exner)*mesh.magSf();
            }}
        }
    }
}
