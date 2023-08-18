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

if (SIgravityWaves || hydrostatic)
{
    // Calculate excess stratification
    volScalarField Gamma = max
    (
        (fvc::grad(thetaf) & kdir)
      - thetaa/(sqr(ocCoeff*runTime.deltaT())*magg),
        dimensionedScalar(theta.dimensions()/dimLength, scalar(0))
    );

    volScalarField G = -sqr(ocCoeff)*runTime.deltaT()*rho*Gamma*
    (
        (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
    );

    thetapf += ocCoeff*runTime.deltaT()*fvc::interpolate(w*Gamma);
    /*
    // Face calculations
    surfaceScalarField Gamma = fvc::snGrad(theta);
    volScalarField G = -sqr(ocCoeff)*runTime.deltaT()*rho*
    (
        kdir & 
        fvc::weightedReconstruct(Gamma*mesh.magSf(), 0.25)
    )*
    (
        kdir & 
        fvc::weightedReconstruct
        (
            gSf/thetaaf + Cp*fvc::snGrad(Exnerp)*mesh.magSf(),
            0.25
        )
    );
    
    thetapf += ocCoeff*runTime.deltaT()*Gamma*phi
                /(rhof*mesh.magSf());
    *///
    if (!Boussinesq) { thetaf = thetaaf + thetapf; }

    // Central matrix coefficient
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (ucorr == 0)
{
    Urhs = -rho*fvc::weightedReconstruct
    (
        offCentre*
        (
            Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()
          + gSf*thetapf/thetaaf
        ),
        0.25
    );
}
solve(uEqn == (Urhs & idir));
solve(wEqn == (Urhs & kdir));

