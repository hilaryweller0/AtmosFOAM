// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  + fvm::div(ImEx*aphi, u, "div(phi,U)")
  + fvc::div((1-ImEx)*aphi, u, "div(phi,U)")
  - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), u)
  //+ 2*rho*(Omega^(U-Ug))
  - (Urhs.oldTime() & idir)
);

// w equation
fvScalarMatrix wEqn
(
    fvm::Sp(rho*muSponge, w)
  - (Urhs.oldTime() & kdir)
);

if (SIgravityWaves || hydrostatic)
{
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);
    volScalarField offCC = fvc::average(offCentre);

    if (!Boussinesq)
    {
        G = max
        (
            -runTime.deltaT()*sqr(offCC)*rho*dThetadz*
            (
                (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
            ),
            dimensionedScalar("", G.dimensions(), scalar(0))
        );
    }
    
    //thetapf += runTime.deltaT()*fvc::interpolate(offCC*w*dThetadz*sign(G));
    //if (!Boussinesq) { thetaf = thetapf + thetaaf; }
    
    // Central matrix coefficient    
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        EulerDdt.fvmDdt(rho, w)
      + fvm::div(ImEx*aphi, w, "div(phi,U)")
      + fvc::div((1-ImEx)*aphi, w, "div(phi,U)")
      - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), w)
    );
}

if (ucorr == 0)
{
    solve(uEqn);
    solve(wEqn);
}
else
{
    solve(uEqn == (Urhs & idir));
    solve(wEqn == (Urhs & kdir));
}

