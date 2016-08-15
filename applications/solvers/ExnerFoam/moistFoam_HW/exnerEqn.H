{
    // Eqn of state and compressibility
    Exner = atmos.ExnerFromTheta(theta);
    Psi == atmos.rho()/Exner;
    rhof = fvc::interpolate(atmos.rho());

    surfaceScalarField G("G", 1/(1+offCentre*dt*muSponge));

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        gradPcoeff = fvc::interpolate(air.Cp()*theta);
        surfaceScalarField gradPcoeff2 = gradPcoeff*G*offCentre*rhof*dt;
        dVdt = -(H.delta() &
        (
            fvc::interpolate(fvc::div(un, rhoU))
         + (2*rhof*(Omega^(Uf-Ug)))
        ));

        // Vu is the flux without g or dpdz (to get the BCs right)
        Vu = G*
        (
            V.oldTime()
          + (1-offCentre)*dt*dVdt.oldTime()
          + offCentre*dt*dVdt
        );

        // add the gravity terms
        V = Vu + G*offCentre*dt*rhof*gd;

        // Flux in the normal direction (including off diagonal part of pressure
        // gradient)

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            U = H.ddirToFlux(V)
              - H.ddirToFluxCorr(gradPcoeff2*fvc::snGrad(Exner)*H.magd());

            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              - (1-offCentre)*fvc::div(U.oldTime())
              + offCentre*fvc::div(U)
              - fvm::laplacian(offCentre*gradPcoeff2, Exner)
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
                U += ExnerEqn.flux()/offCentre;
                V -= gradPcoeff2*fvc::snGrad(Exner)*H.magd();
                rho = Exner*Psi;
                rhof = fvc::interpolate(rho);
                rhoU = H.reconstructd(V);
                un = U/rhof;
                divu = fvc::div(un);
                Uf = fvc::interpolate(rhoU/rho);
                dVdt += rhof*gd
                     - H.magd()*rhof*gradPcoeff*fvc::snGrad(Exner)
                     - muSponge*V;
            }
        }
    }
}
