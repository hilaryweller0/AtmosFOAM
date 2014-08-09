{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    surfaceScalarField G("G", 1+offCentre*dt*muSponge);

    if (SIgravityWaves)
    {
        G =1/(G - sqr(dt*offCentre)*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner));
        vperp = Uf - (Uf & H.delta())*H.delta()/sqr(H.magd());
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
        dVdt = -H.delta() &
        (
            fvc::interpolate(fvc::div(U,u)) 
         + (2*rhof*(Omega^(Uf-Ug)))
        );
        nonlinearMomentum = fvc::div(U,u);

        // Vu is the flux without g or dpdz (to get the BCs right)
        Vu = G*(V.oldTime() + offCentre*dt*dVdt);

        // add the gravity terms
        V = Vu + G*offCentre*dt*rhof*gd;

        // Flux in the normal direction (including off diagonal part of pressure
        // gradient)
        surfaceScalarField gradPcoeff = G*offCentre*Cp*dt
                                        *rhof*thetaf*H.magd();

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            U = H.ddirToFlux(V)
              - H.ddirToFluxOffDiag(gradPcoeff*fvc::snGrad(Exner));

            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, Exner)
              + (1-offCentre)*divU
              + offCentre*fvc::div(U)
              - fvm::laplacian
                (
                    offCentre*H.Hdiag()*gradPcoeff/mesh.magSf(),
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
                U += ExnerEqn.flux()/offCentre;
                V -= gradPcoeff*fvc::snGrad(Exner);
                u = H.reconstructd(V/rhof);
                Uf = fvc::interpolate(u);
            }
        }
    }
    divUf = fvc::div(Uf & mesh.Sf());
}

