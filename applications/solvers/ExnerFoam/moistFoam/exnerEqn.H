{
    kappam = Rm/cpml;
    
    Psi = rho*pow(rho*R*thetaRho/pRef, kappam/(kappam-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);
    
    surfaceScalarField G("G", 1/(1+offCentre*dt*muSponge));

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
        gradPcoeff = R*fvc::interpolate(thetaRho/kappam);
        gradPcoeff2 = gradPcoeff*G*offCentre*rhof*dt;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            U = H.ddirToFlux(V)
              - H.ddirToFluxCorr(gradPcoeff2*fvc::snGrad(Exner)*H.magd());

            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              - (1-offCentre)*drhodt.oldTime()
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
                //Uv = U/rhof;
                V -= gradPcoeff2*fvc::snGrad(Exner)*H.magd();
                u = H.reconstructd(V/rhof);
                Uf = fvc::interpolate(u);
            }
        }
    }
}
