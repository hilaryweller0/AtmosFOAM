{
    kappam = Rm/cpml;
    
    // Eqn of state and compressibility
    Exner = moistGasExner(rho, thetaRho, kappam, R, pRef);
    Psi == rho/Exner;
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

        // add the gravity terms and pressure gradient terms
        gradPcoeff = gradPCoeff(cpml, thetaRho, rv, epsilon);
        V = Vu + G*offCentre*dt*rhof*
            (gd - gradPcoeff*fvc::snGrad(Exner)*H.magd());

        U = H.ddirToFlux(V);
        u = H.reconstructd(V/rhof);
        Uf = fvc::interpolate(u);
    }
}

