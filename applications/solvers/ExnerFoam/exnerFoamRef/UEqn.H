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
    // Special Imex for wEqn when SIgravityWaves
    surfaceScalarField ImExw = ImEx;
    if (SIgravityWaves)
    {
        volScalarField Ndt = runTime.deltaT()*Foam::sqrt(max
        (
            -(g & fvc::grad(thetaf))/theta,
            dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
        ));

        ImExw += (1 - ImExw)*0.5*(sign(maxInterp.interpolate(Ndt) - 1)+1);
    }

    wEqn += fvScalarMatrix
    (
        EulerDdt.fvmDdt(rho, w)
      + fvm::div(ImExw*aphi, w, "div(phi,U)")
      + fvc::div((1-ImExw)*aphi, w, "div(phi,U)")
      - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), w)
    );
}

if (SIgravityWaves || hydrostatic)
{
    volScalarField dThetadz = max
    (
        fvc::grad(theta) & kdir,
        dimensionedScalar("", theta.dimensions()/dimLength, scalar(0))
    );
    volScalarField offCC = fvc::average(offCentre);

    volScalarField G = -runTime.deltaT()*sqr(offCC)*rho*dThetadz*
    (
        (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
    );
    
    // Central matrix coefficient
    wEqn += fvScalarMatrix(fvm::Sp(G, w));

    // Temporary storage of thetaf increment for after solution of wEqn
    thetaf = runTime.deltaT()*fvc::interpolate(offCC*w*dThetadz);
    //thetapf += runTime.deltaT()*fvc::interpolate(offCC*w*dThetadz);
    //if (!Boussinesq) { thetaf = thetaaf + thetapf; }
    
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

if (SIgravityWaves || hydrostatic)
{
    // thetepf and thetaf incremented based on temporary increment stored in
    // thetaf
    thetapf += thetaf;
    thetaf = thetaaf + thetapf;
}

