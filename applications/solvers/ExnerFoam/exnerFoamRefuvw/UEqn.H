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
    
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);
    
    volScalarField G
    (
        "G",
        -adt*rho*dThetadz*
        (
            (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
        )
    );
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    surfaceScalarField Gf = linearInterpolate(G);

    //surfaceScalarField Gf = -adt*rhof*fvc::snGrad(theta)
    //                    *(Cp*fvc::snGrad(Exnerp) + gSf/(mesh.magSf()*thetaaf));

    // Apply this only where on vertical pointing faces where G>0
    surfaceScalarField applySIg = mag(gUnitNormal)*max(sign(Gf), scalar(0));
    
    //thetapf +=  applySIg*adt*phi*fvc::snGrad(theta)/(rhof*mesh.magSf());
    thetapf +=  applySIg*adt*fvc::interpolate(w*dThetadz);
    thetaf = thetapf + thetaaf;

    //volScalarField G("G", fvc::average(Gf));
    //G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    
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
    volVectorField Urhs = -rho*
        fvc::reconstruct(Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf());
    
    solve(uEqn == (Urhs & idir));
    
    Urhs -= rho*fvc::reconstruct(thetapf*gSf/thetaaf);

    solve(wEqn == (Urhs & kdir));
}

