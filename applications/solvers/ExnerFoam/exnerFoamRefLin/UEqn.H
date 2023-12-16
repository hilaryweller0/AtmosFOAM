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
    EulerDdt.fvmDdt(rho, w)
  + fvm::div(ImEx*aphi, w, "div(phi,U)")
  + fvc::div((1-ImEx)*aphi, w, "div(phi,U)")
  - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), w)
  + fvm::Sp(rho*muSponge, w)
  - (Urhs.oldTime() & kdir)
);

if (SIgravityWaves)
{
    volScalarField G = runTime.deltaT()*rho*N2*fvc::average(sqr(offCentrePG));

    // Central matrix coefficient
    wEqn += fvScalarMatrix(fvm::Sp(G, w));

    thetapf += offCentrePG*runTime.deltaT()*fvc::interpolate(w*dThetaadz);
}

if (ucorr == 0)
{
    Urhs = -rho*fvc::weightedReconstruct
    (
        offCentrePG*
        (
            Cp*thetaaf*fvc::snGrad(Exnerp)*mesh.magSf()
          + gSf*thetapf/thetaaf
        ),
        0.25
    );
}
solve(uEqn == (Urhs & idir));
solve(wEqn == (Urhs & kdir));

