// Implicit solution of the Momentum equation for the turbulence

// Separate u and w equations
fvScalarMatrix uEqn
(
    fvm::ddt(rho, u)
  + fvm::div(phi, u, "div(phi,U)")
  - fvm::laplacian(rho*turbulence->nuEff(), u)
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
    volScalarField G("G", fvc::average(-adt*rhof*Cp*gradnTheta*fvc::snGrad(Exner)));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetaf += min(adt*wf, 1/mesh.deltaCoeffs())*gradnTheta;

    /*// Cell center calculations
    volScalarField dThetadz("dThetadz", fvc::grad(theta) & kdir);
    volScalarField G(-fvc::average(adt)*rho*Cp*dThetadz*(fvc::grad(Exner) & kdir));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    thetaf += adt*fvc::interpolate(w*dThetadz);*/

    // Central matrix coefficient    
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        fvm::ddt(rho, w)
      + fvm::div(phi, w, "div(phi,U)")
      - fvm::laplacian(rho*turbulence->nuEff(), w)
    );
}

if (implicitU)
{
    volVectorField Urhs = rho*fvc::weightedReconstruct
    (
        gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf(), 0.25
    );

    solve(uEqn == (Urhs & idir));
    solve(wEqn == (Urhs & kdir));
}

