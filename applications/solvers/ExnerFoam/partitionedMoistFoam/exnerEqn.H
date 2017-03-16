{
    // Global Eqn of state and from each partition
    Exner = atmosParts.exnerFromState();
    Exner.write();
    for(label ip = 0; ip < atmosParts.size(); ip++)
    {
        volScalarField Exneri
        (
            atmosParts[ip].partitionName()+"Exner",
            atmosParts[ip].exnerFromState()
        );
        Exneri.write();
    }
    FatalErrorIn("partitionedMoistFoam") << exit(FatalError);
/*
    // Implicit parts of momentum equation
    surfaceScalarField G("G", 1/(1+offCentre*dt*muSponge));
    
    for (int corr=0; corr<nCorr; corr++)
    {
        // Compressibility and updated density
        if (corr != 0)
        {
            #include "rhoEqn.H"
        }
        Psi == rho/Exner;
        rhof = fvc::interpolate(rho);

        // Explicit part of rate of change of flux
        gradPcoeff = fvc::interpolate
        (
            air.Cp()*theta*atmos.rhoR()
           /(atmos.sumDensity()*air.R()*atmos.volGas())
        );
        surfaceScalarField gradPcoeff2 = gradPcoeff*G*offCentre*rhof*dt;
        dVdt = -(H.delta() &
        (
            fvc::interpolate(fvc::div(U, u))
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
              - (1-offCentre)*dRhodt.oldTime()
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
                u = H.reconstructd(V/rhof);
                divu = fvc::div(U/rhof);
                Uf = fvc::interpolate(u);
                dVdt += rhof*gd
                     - H.magd()*rhof*gradPcoeff*fvc::snGrad(Exner)
                     - muSponge*V;
            }
        }
    }*/
}

