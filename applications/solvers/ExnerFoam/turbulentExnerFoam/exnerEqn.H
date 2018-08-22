{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dPhidt = -mesh.Sf() & fvc::interpolate
        (
            fvc::div(phi, U)
          + ((turbulence->divDevRhoReff(U)) & U)
        );

        // phiU is the flux without g or dpdz (to get the BCs right)
        phiU = phi.oldTime() + dt/(1+ocCoeff)*
        (
             dPhidt
           + ocCoeff*dPhidt.oldTime()
        );

        // add the gravity terms (and the non-orthogonal grad Exner terms)
        phi = phiU + dt/(1+ocCoeff)*rhof*gSf;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + fvc::div(phi)
              - fvm::laplacian(dt/(1+ocCoeff)*Cp*rhof*thetaf, Exner)
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
                phi -= dt/(1+ocCoeff)*Cp*rhof*thetaf*fvc::snGrad(Exner)
                     *mesh.magSf();
                U = fvc::reconstruct(phi/rhof);
            }
        }
        // Update the pressure and temperature based on the new Exner
        thermo.p() = pRef*pow(Exner, 1/kappa);
        thermo.T() = theta*Exner;
        thermo.he() == thermo.he(thermo.p(),thermo.T());
        thermo.correct();
    }
}

