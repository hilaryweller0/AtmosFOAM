{
    rho = pRef/(R*theta)*pow(Exner, (1-kappa)/kappa);
    surfaceScalarField rhof = fvc::interpolate(rho);
    for (int corr=0; corr<nCorr; corr++)
    {
        volScalarField Psi = (1-kappa)/kappa*rho/Exner;

        // Rate of change of velocity
        dVolFluxDt = 
        (
            mesh.Sf() &
            (
              -fvc::interpolate
               (
                   fvc::div(volFlux, u, "u") - u * fvc::div(volFlux)
                   - ((turbulence->divDevRhoReff(u)) & u)/rho,
                   "convection"
               )
            )
        )
        + gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf();
        
        // Coriolis term
        if (magSqr(Omega.value()) > VSMALL)
        {
            dVolFluxDt -= mesh.Sf() & (2*(Omega ^ Uf));
        }
        
        // Velocity predictor
        volFlux = volFlux.oldTime()
                + dt*(1-offCentre)*dVolFluxDt.oldTime()
                + dt*offCentre*dVolFluxDt;
        
        flux = rhof*volFlux;

        // Helmholtz equation
        for (int nonOrth=0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            // First guess at ExnerPrime
            ExnerPrime *= 0;
            fvScalarMatrix ExnerEqn
            (
                fvm::Sp(Psi/dt, ExnerPrime)
              + fvc::ddt(rho)
              + (1-offCentre)*fvc::div(flux.oldTime())
              + offCentre*fvc::div(flux)
              - fvm::laplacian(dt*sqr(offCentre)*Cp*rhof*thetaf, ExnerPrime)
            );

            ExnerEqn.solve(Exner.name());
            Info << "Time " << runTime.timeName() << " ExnerPrime convergence "
                 << (fvc::domainIntegrate(mag(ExnerPrime))/sum(mesh.V())).value()
                 << endl;
        }

        // Back-substitutions
        Exner += ExnerPrime;
        dVolFluxDt -= Cp*thetaf*fvc::snGrad(ExnerPrime)*mesh.magSf();
        volFlux -= dt*offCentre*Cp*thetaf*fvc::snGrad(ExnerPrime)*mesh.magSf();
        
        rho += Psi*ExnerPrime;
//        rho  = pRef/(R*theta)*pow(Exner, (1-kappa)/kappa);
        rhof = fvc::interpolate(rho);
        flux = rhof*volFlux;
        u = fvc::reconstruct(volFlux);
        Uf = fvc::interpolate(u, "Uf");
        
        // Update the pressure and temperature based on the new Exner
        thermo.p() = pRef*pow(Exner, 1/kappa);
        thermo.T() = theta*Exner;
        thermo.he() == thermo.he(thermo.p(),thermo.T());
        thermo.correct();
    }
}

