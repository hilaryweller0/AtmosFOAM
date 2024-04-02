// Implicit solution of the Momentum equation for the turbulence

// Separate u, v and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  + fvm::div(ImEx*aphi, u, "div(phi,U)")
  + fvc::div((1-ImEx)*aphi, u, "div(phi,U)")
  - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), u)
  //+ 2*rho*(Omega^(U-Ug))
  - (Urhs.oldTime() & idir)
);

fvScalarMatrix vEqn
(
    EulerDdt.fvmDdt(rho, v)
  + fvm::div(ImEx*aphi, v, "div(phi,U)")
  + fvc::div((1-ImEx)*aphi, v, "div(phi,U)")
  - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), v)
  //+ 2*rho*(Omega^(U-Ug))
  - (Urhs.oldTime() & jdir)
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
    /*volScalarField dThetadz = max
    (
        (fvc::grad(thetaf) & kdir),
        dimensionedScalar(theta.dimensions()/dimLength, scalar(0))
    );*/

    volScalarField G = -runTime.deltaT()*rho*dThetaadz*
    (
        fvc::weightedReconstruct
        (
            offCentrePG*offCentreGW
                *(gSf/thetaaf/* + Cp*fvc::snGrad(Exnerp)*mesh.magSf()*/),
            0.25
        )
        & kdir
    );

    // Central matrix coefficient
    wEqn += fvm::Sp(G, w);

    thetapf += offCentreGW*runTime.deltaT()*fvc::interpolate(w*dThetaadz);
    //thetaf = thetaaf + thetapf;
}

if (ucorr == 0)
{
    Urhs = -rho*fvc::weightedReconstruct
    (
        offCentrePG*
        (
            Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()
          + gSf*thetapf/thetaaf
        ),
        0.25
    );
}
if (mesh.solutionD()[0] == 1) solve(uEqn == (Urhs & idir));
if (mesh.solutionD()[1] == 1) solve(vEqn == (Urhs & jdir));
if (mesh.solutionD()[2] == 1) solve(wEqn == (Urhs & kdir));

