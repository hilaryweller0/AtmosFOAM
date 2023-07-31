// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  + fvm::div(phi, u, "div(phi,U)")
  - fvm::laplacian(ocAlpha*rhof*fvc::interpolate(turbulence->nuEff()), u)
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
    surfaceScalarField thetaInc = -ocAlpha*runTime.deltaT()*fvc::snGrad(theta)*gUnitNormal;
    thetaInc = max(thetaInc, dimensionedScalar(thetaInc.dimensions(), scalar(0)));

    volScalarField G = -fvc::average(thetaInc)*rho*Cp*(fvc::grad(Exnerg) & kdir);

    thetaf +=  fvc::interpolate(w)*thetaInc;

    // Central matrix coefficient    
    wEqn += fvScalarMatrix(fvm::Sp(ocCoeff*G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        EulerDdt.fvmDdt(rho, w)
      + fvm::div(phi, w, "div(phi,U)")
      - fvm::laplacian(ocAlpha*rho*turbulence->nuEff(), w)
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

