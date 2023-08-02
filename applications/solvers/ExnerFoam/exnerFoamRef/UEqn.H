// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), u)
  //+ 2*rho*(Omega^(U-Ug))
  - (Urhs.oldTime() & idir)
);
if (advectionType == advType::Implicit)
{
    uEqn += fvScalarMatrix(fvm::div(aphi, u, "div(phi,U)"));
}
else
{
    uEqn += fvc::div(aphi, u, "div(phi,U)");
}

// w equation
fvScalarMatrix wEqn
(
    fvm::Sp(rho*muSponge, w)
  - (Urhs.oldTime() & kdir)
);

if (SIgravityWaves || hydrostatic)
{
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);

    if (!Boussinesq)
    {
        volScalarField G
        (
            -runTime.deltaT()*fvc::average(offCentre)*rho*dThetadz*
            (
                (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
            )
        );
        G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    }
    
    thetapf += offCentre*runTime.deltaT()*fvc::interpolate(w*dThetadz)
                *sign(fvc::interpolate(G));

    if (!Boussinesq)
    {
        thetaf = thetapf + thetaaf;
    }
    
    // Central matrix coefficient    
    wEqn += fvScalarMatrix(fvm::Sp(ocCoeff*G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        EulerDdt.fvmDdt(rho, w)
      - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), w)
    );
    if (advectionType == advType::Implicit)
    {
        wEqn += fvScalarMatrix(fvm::div(aphi, w, "div(phi,U)"));
    }
    else
    {
        wEqn += fvc::div(aphi, w, "div(phi,U)");
    }
}

if (advectionType == advType::Implicit)
{
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
}

