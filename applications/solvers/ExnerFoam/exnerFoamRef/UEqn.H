// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    EulerDdt.fvmDdt(rho, u)
  + fvc::div(phiMid, u, "div(phi,U)")
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

if (SIgravityWaves || hydrostatic)
{
    if (!Boussinesq)
    {
        gradnTheta = fvc::snGrad(theta);
    
        // Face calculations
        G = fvc::average
        (
            -offCentre*runTime.deltaT()*rhof*gradnTheta*
            (
                Cp*fvc::snGrad(Exnerp) + gSf/(mesh.magSf()*thetaaf)
            )
        );
        G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    }
    
    thetapf += min
    (
        runTime.deltaT()*offCentre*phi/(rhof*mesh.magSf()),
        1/mesh.deltaCoeffs()
    )*gradnTheta;

    /*// Cell center calculations
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);

    volScalarField G
    (
        -fvc::average(adt)*rho*dThetadz*
        (
            (g&kdir)/thetaa + (Cp*fvc::grad(Exnerp) & kdir)
        )
    );
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetapf += adt*fvc::interpolate(w*dThetadz);*/

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
      + fvc::div(phiMid, w, "div(phi,U)")
      - fvm::laplacian(offCentre*rhof*fvc::interpolate(turbulence->nuEff()), w)
    );
}

//if (advectionType == advType::Implicit)
//{
    if (ucorr == 0)
    {
        Urhs = rho*fvc::weightedReconstruct
        (
            offCentre
           *(Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf() - gSf*thetapf/thetaaf),
           0.25
        );
    }
    
    solve(uEqn == (Urhs & idir), "UEqn");
    solve(wEqn == (Urhs & kdir), "UEqn");
//}

