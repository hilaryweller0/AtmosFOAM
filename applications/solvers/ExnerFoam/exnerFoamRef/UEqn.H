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
    surfaceScalarField gradnTheta(fvc::snGrad(theta));

    surfaceScalarField Ndt
    (
        Foam::sqrt
        (
            max(-sqr(runTime.deltaT())*gSfhat*gradnTheta/thetaf, VSMALL)
        )
    );
    surfaceScalarField wf(phi/(rhof*mesh.magSf()));

    surfaceScalarField adt = runTime.deltaT()*max(scalar(0), 1-1/Ndt);
    
    // Face calculations
    volScalarField G
    (
        "G",
        fvc::average
        (
            -adt*rhof*gradnTheta*
            (
                Cp*fvc::snGrad(Exnerp) + gSf/(mesh.magSf()*thetaaf)
            )
        )
    );
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetapf += min(adt*wf, 1/mesh.deltaCoeffs())*gradnTheta;

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

    thetaf = thetapf + thetaaf;
    // Central matrix coefficient    
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
    volVectorField Urhs = rho*fvc::weightedReconstruct
    (
        Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()
      + gSf*thetapf/thetaaf,
        1
    );

    solve(uEqn == (Urhs & idir));
    solve(wEqn == (Urhs & kdir));
}

