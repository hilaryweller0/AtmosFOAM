// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  + fvm::div(phiNew, u, "div(phi,u)")
  + fvc::div(phiOld, u.oldTime(), "div(phi,u)")
  //- fvm::laplacian(rho*turbulence->nuEff(), u)
  //+ 2*rho*(Omega^(U-Ug))
  - (1-ocAlpha)*(Urhs.oldTime() & idir)
);

// w equation
fvScalarMatrix wEqn
(
    fvm::Sp(rho*muSponge, w)
  - (1-ocAlpha)*(Urhs.oldTime() & kdir)
);

if (SIgravityWaves || hydrostatic)
{
    surfaceScalarField gradnTheta(fvc::snGrad(theta));

    surfaceScalarField Ndt
    (
        Foam::sqrt
        (
            max(-sqr(runTime.deltaT())*gSfhat*gradnTheta/thetaf, VSMALL)
        )
    );
    surfaceScalarField wf(phi/(rhof*mesh.magSf()));

    surfaceScalarField adt = runTime.deltaT()*max(ocAlpha, 1-1/Ndt);
    
    // Face calculations
    volScalarField G("G", fvc::average(-adt*rhof*Cp*gradnTheta*fvc::snGrad(Exner)));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetaf += min(adt*wf, 1/mesh.deltaCoeffs())*gradnTheta;

    /*// Cell center calculations
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);
    volScalarField G(-fvc::average(adt)*rho*Cp*dThetadz*(fvc::grad(Exner) & kdir));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetaf += adt*fvc::interpolate(w*dThetadz);*/

    // Central matrix coefficient    
    wEqn += fvScalarMatrix(fvm::Sp(ocCoeff*G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        EulerDdt.fvmDdt(rho, w)
      + fvm::div(phiNew, w, "div(phi,w)")
      + fvc::div(phiOld, w.oldTime(), "div(phi,w)")
      //- fvm::laplacian(rho*turbulence->nuEff(), w)
    );
}

if (implicitU)
{
    solve(uEqn == ocAlpha*(Urhs & idir));
    solve(wEqn == ocAlpha*(Urhs & kdir));
}

