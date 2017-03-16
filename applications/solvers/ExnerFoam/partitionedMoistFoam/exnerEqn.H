{
    // Global Eqn of state and from each partition
    Exner = atmosParts.exnerFromState();

    // Implicit parts of momentum equation
    //surfaceScalarField G("G", 1/(1+offCentre*dt*muSponge));
    
    for (int corr=0; corr<nCorr; corr++)
    {
        volScalarField Psi = atmosParts.rho()/Exner;

        // Contributions from each partition
        
        // Diffusivity in each partition
        surfaceScalarField gradPcoeff
        (
            IOobject("gradPcoeff", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("d",dimensionSet(1,-1,-2,0,0), scalar(0))
        );
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            gradPcoeff += air.Cp()/air.R()*fvc::interpolate
            (
                parti.sigmaRho()*parti.sigma()*parti.rhoR()*parti.theta()
              / (parti.sigmaRho()*parti.volGas())
            );
        }

/*        surfaceScalarField gradPcoeff2 = gradPcoeff*G*offCentre*rhof*dt;
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
        }*/
    }
}

