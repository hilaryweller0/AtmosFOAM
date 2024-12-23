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
    // Calculate stratification
    volScalarField Gamma = max
    (
        (fvc::grad(thetaf) & kdir),
        dimensionedScalar(theta.dimensions()/dimLength, scalar(0))
    );

    volScalarField G = -runTime.deltaT()*rho*Gamma*
    (
        fvc::weightedReconstruct
        (
            offCentrePG*offCentreGW*Cp*fvc::snGrad(Exner)*mesh.magSf(),
            0.25
        )
        & kdir
    );

    // Central matrix coefficient
    wEqn += fvScalarMatrix(fvm::Sp(G, w));

    thetaf += offCentreGW*runTime.deltaT()*fvc::interpolate(w*Gamma);
}

if (ucorr == 0)
{
    Urhs = rho*fvc::weightedReconstruct
    (
        offCentrePG*(gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()), 0.25
    );
}
solve(uEqn == (Urhs & idir));
solve(wEqn == (Urhs & kdir));

