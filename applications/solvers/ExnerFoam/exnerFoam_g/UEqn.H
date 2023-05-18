// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    fvm::ddt(rho, u)
  + fvm::div(phi, u)
  //- fvm::laplacian(rho*turbulence->nuEff(), u)
  //+ 2*rho*(Omega^(U-Ug))
);

// w equation
fvScalarMatrix wEqn(fvm::Sp(rho*muSponge, w));

if (SIgravityWaves || hydrostatic)
{
    const dimensionedScalar adt = ocAlpha*runTime.deltaT();
    
    surfaceScalarField Gf = -adt*rhof*Cp*gradTheta*fvc::snGrad(Exnerg);

    // Apply this only where on vertical pointing faces where G>0
    surfaceScalarField applySIg = mag(gUnitNormal)*max(sign(Gf), scalar(0));
    
    thetaf +=  applySIg*adt*gradTheta*phi/(rhof*mesh.magSf());

    volScalarField G("G", fvc::average(Gf));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        fvm::ddt(rho, w)
      + fvm::div(phi, w)
      //- fvm::laplacian(rho*turbulence->nuEff(), w)
    );
}

if (implicitU)
{
    volVectorField Urhs = rho*fvc::reconstruct
    (
        (ghGradThetaByTheta - Cp*thetaf*fvc::snGrad(Exnerg))*mesh.magSf()
    );
    solve(uEqn == (Urhs & idir));
    solve(wEqn == (Urhs & kdir));
}

