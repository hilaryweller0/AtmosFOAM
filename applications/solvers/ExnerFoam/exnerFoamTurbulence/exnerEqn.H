{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    thetaf = fvc::interpolate(theta);

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dPhidt = -mesh.Sf() & fvc::interpolate
        (
//              (UEqn & U) - fvc::ddt(rho, U)
            fvc::div(phi, U)
          + ((turbulence->divDevRhoReff(U)) & U)
          + 2*rho*(Omega^(U-Ug))
          - muSponge*rho*(sqr(ghat) & U)
        );

        // phiU is the flux without g or dpdz (to get the BCs right)
        phiU = phi.oldTime()
             + (1-offCentre)*dt*dPhidt.oldTime()
             + offCentre*dt*dPhidt;

        // add the gravity terms (and the non-orthogonal grad Exner terms)
        phi = phiU + offCentre*dt*rhof*gSf;

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
                ExnerEqn.solve(mesh.solver(Exner.name() + "Final"));
            }
            else
            {
                ExnerEqn.solve(mesh.solver(Exner.name()));
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi -= offCentre*dt*Cp*rhof*thetaf*fvc::snGrad(Exner)*mesh.magSf();
                //phi += ExnerEqn.flux();
                //U = H.reconstruct(phi/rhof);
                U = fvc::reconstruct(phi/rhof);
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

