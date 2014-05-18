{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Psi*Exner;
    rhof = fvc::interpolate(rho);

    surfaceScalarField G("G", 1+offCentre*dt*muSponge);

    if (SIgravityWaves)
    {
        G =1/(G - sqr(dt*offCentre)*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner));
        if (ucorr == 0 && corr == 0)
        //if (ucorr < nOuterCorr-1 || corr < nCorr-1)
        {
            vperp = rhof*(Uf - (Uf & H.delta())*H.delta()/sqr(H.magd()));
            thetaf = thetaf.oldTime()
               - offCentre*dt*(vperp & fvc::interpolate(fvc::grad(theta)))/rhof;
        }
        else
        {
            thetaf = fvc::interpolate(theta)
                   + offCentre*dt*V/(rhof*H.magd())*fvc::snGrad(theta);
        }
    }
    else
    {
        thetaf = fvc::interpolate(theta);
        G = 1/G;
    }

    // Explicit part of rate of change of flux
    dVdt = -H.delta() & (2*rhof*(Omega^(Uf-Ug)));

    // Vu is the flux without g or dpdz (to get the BCs right)
    Vu = G*(V.oldTime() + offCentre*dt*dVdt);

    // add the gravity terms
    V = Vu + G*offCentre*dt*rhof*gd;

    // Flux in the normal direction (including off diagonal part of pressure
    // gradient)
    surfaceScalarField gradPcoeff = G*offCentre*Cp*dt
                                    *rhof*thetaf*H.magd();

    U = H.ddirToFlux(V)
      - H.ddirToFluxOffDiag(gradPcoeff*fvc::snGrad(Exner));

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            fvm::ddt(Psi, Exner)
          + (1-offCentre)*divU
          + offCentre*fvc::div(U)
          - fvm::laplacian(offCentre*H.Hdiag()*gradPcoeff/mesh.magSf(), Exner)
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

//{
//    surfaceVectorField vecg
//    (
//        "vecg", fvc::interpolate(H.reconstruct(H.ddirToFlux(gd)))
//    );
//    vecg += (H.ddirToFlux(gd) - (vecg & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
//    vecg.write();
//    
//    thetaf.write();
//    surfaceVectorField gradExner
//    (
//        "gradExner",
//        fvc::interpolate(H.reconstruct(H.ddirToFlux(thetaf*fvc::snGrad(Exner))))
//    );
//    gradExner += (H.ddirToFlux(thetaf*fvc::snGrad(Exner)) - (gradExner & mesh.Sf()))
//                 *mesh.Sf()/sqr(mesh.magSf());
//    gradExner.write();
//}

