{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    surfaceScalarField G("G", 1+offCentre*dt*muSponge);

    if (SIgravityWaves)
    {
        G =1/(G - sqr(dt*offCentre)*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner));
        vperp = Uf - (Uf & mesh.Sf())*mesh.Sf()/sqr(mesh.magSf());
        thetaf = thetaf.oldTime()
               - offCentre*dt*(vperp & fvc::interpolate(fvc::grad(theta)));
    }
    else
    {
        thetaf = fvc::interpolate(theta);
        G = 1/G;
    }

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dPhidt = -mesh.Sf() &
        (
            fvc::interpolate(fvc::div(phi,U)) 
         + (2*rhof*(Omega^(Uf-Ug)))
        );

        // phiU is the flux without g or dpdz (to get the BCs right)
        phiU = G*(phi.oldTime() + offCentre*dt*dPhidt);

        // add the gravity terms (and the non-orthogonal grad Exner terms)
        phi = phiU + G*offCentre*dt*rhof*gSf;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + (1-offCentre)*divPhi
              + offCentre*fvc::div(phi)
              - fvm::laplacian
                (
                    G*sqr(offCentre)*dt*Cp*rhof*thetaf,
                    Exner
                )
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
                phi += ExnerEqn.flux()/offCentre;
                //U = H.reconstruct(phi/rhof);
                U = fvc::reconstruct(phi/rhof);
                Uf = fvc::interpolate(U);
            }
        }
    }
}

