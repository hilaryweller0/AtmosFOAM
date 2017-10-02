{
    Psi == pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    //rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    surfaceScalarField G("G", 1+offCentre*dt*muSponge);

    if (SIgravityWaves)
    {
        G ==1/(G - sqr(dt*offCentre)*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner));
        thetaf = thetaf.oldTime()
               - offCentre*dt*
                 (
                    (Uf & fvc::interpolate(fvc::grad(theta)))
                 - 0.5*V/(rhof*H.magd())*fvc::snGrad(theta)
                 );
    }
    else
    {
        thetaf = fvc::interpolate(theta);
        G == 1/G;
    }

    surfaceScalarField gradPcoeff = dt*Cp*rhof*thetaf*G*offCentre;

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dVdt = -H.delta() &
        (
            fvc::interpolate(fvc::div(U,u))
         + (2*rhof*(Omega^(Uf-Ug)))
        );

        // Vu is the flux without g or dpdz (to get the BCs right)
        Vu = G*(V.oldTime() + offCentre*dt*dVdt);

        // add the gravity terms
        V = Vu + G*offCentre*dt*rhof*gd;

        // Flux in the normal direction (including off diagonal part of pressure
        // gradient)
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            U = H.ddirToFlux(V)
              - H.ddirToFluxCorr(gradPcoeff*fvc::snGrad(Exner)*H.magd());

            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + (1-offCentre)*divU
              + offCentre*fvc::div(U)
              - fvm::laplacian(offCentre*gradPcoeff, Exner)
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
                V -= gradPcoeff*fvc::snGrad(Exner)*H.magd();
                u = H.reconstructd(V)/rho;
                Uf = fvc::interpolate(u);
            }
        }
    }
}
