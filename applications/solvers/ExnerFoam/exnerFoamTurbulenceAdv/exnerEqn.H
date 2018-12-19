{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    thetaf = fvc::interpolate(theta);

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dvolFluxdt = -mesh.Sf() & fvc::interpolate
        (
            fvc::div(volFlux, U)
          - U*fvc::div(volFlux)
//          + ((turbulence->divDevRhoReff(U)) & U)/rho
          + 2*(Omega^(U-Ug))
          - muSponge*(sqr(ghat) & U)
        );

        // volFluxU is the flux without g or dpdz (to get the BCs right)
        volFluxU = volFlux.oldTime()
             + (1-offCentre)*dt*dvolFluxdt.oldTime()
             + offCentre*dt*dvolFluxdt;

        // add the gravity terms
        volFlux = volFluxU + offCentre*dt*gSf;
        phi = volFlux*rhof;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + fvc::div(phi)
              - fvm::laplacian(Cp*rhof*thetaf*offCentre*dt, Exner)
            );

            if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
            {
                ExnerEqn.solve(Exner.name() + "Final");
            }
            else
            {
                ExnerEqn.solve(Exner.name());
            }

            if (nonOrth == nNonOrthCorr)
            {
                volFlux -= offCentre*dt*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf();
                
                phi = volFlux*rhof;
                U = fvc::reconstruct(volFlux);
                Uf = fvc::interpolate(U);
            }
        }
        // Update the pressure and temperature based on the new Exner
        thermo.p() = pRef*pow(Exner, 1/kappa);
        thermo.T() = theta*Exner;
        thermo.he() == thermo.he(thermo.p(),thermo.T());
        thermo.correct();
    }
}

