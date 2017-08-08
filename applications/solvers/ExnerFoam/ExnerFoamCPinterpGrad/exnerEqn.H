{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    thetaf = thetaf.oldTime() - dt*(Uf & fvc::interpolate(fvc::grad(thetaf)));
    bf = thetaf * g.unitFaceNormal();
    b = fvc::reconstruct(bf * mesh.magSf());
    theta == (b & g.unit());
    thetaf = mag(g.unitFaceNormal())*thetaf + (1.0 - mag(g.unitFaceNormal()))*fvc::interpolate(theta, "thetaFromb");

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dVdt = -H.delta() &
        (
            fvc::interpolate(fvc::div(U,u))
         + (2*rhof*(Omega^(Uf-Ug)))
        );

        // Vu is the flux without g or dpdz (to get the BCs right)
        Vu = G*(V.oldTime() + dt*dVdt);

        // add the gravity terms
        V = Vu + G*dt*rhof*gd;

        // Flux in the normal direction (including off diagonal part of pressure
        // gradient)
        surfaceScalarField gradPcoeff("gradPcoeff", G*Cp*dt*rhof*thetaf);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            U = H.ddirToFlux(V)
              - H.ddirToFluxCorr(gradPcoeff*fvc::snGrad(Exner)*H.magd());

            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + fvc::div(U)
              - fvm::laplacian(gradPcoeff, Exner)
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
                U += ExnerEqn.flux();
                V -= gradPcoeff*fvc::snGrad(Exner)*H.magd();
                u = H.reconstructd(V)/rho;
                Uf = fvc::interpolate(u);
            }
        }
    }
}
