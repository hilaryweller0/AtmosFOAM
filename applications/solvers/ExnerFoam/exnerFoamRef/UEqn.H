// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + fvm::laplacian(rho*turbulence->nuEff(), U)
  + 2*rho*(Omega^(U-Ug))
);

/*// w equation
fvScalarMatrix wEqn(fvm::Sp(rho*muSponge, w));
*/

/*if (SIgravityWaves || hydrostatic)
{
    dimensionedScalar adt = ocAlpha*runTime.deltaT();

    surfaceScalarField Gf = -adt*rhof*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner);

    // Apply this only where on vertical pointing faces where G>0
    surfaceScalarField applySIg = mag(gUnitNormal)*max(sign(Gf), scalar(0));
    
    thetaf += applySIg*adt*phi*fvc::snGrad(theta)/(rhof*mesh.magSf());

    volScalarField G = fvc::average(Gf);
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}*/

/*if (!hydrostatic)
{
    wEqn += fvScalarMatrix
    (
        fvm::ddt(rho, w)
      + fvm::div(phi, w)
      - fvm::laplacian(rho*turbulence->nuEff(), w)
    );
}
*/
if (implicitU)
{
    solve
    (
        UEqn == rho*
        (
            g*thetap/thetaa
          - Cp*(thetaa + thetap)*fvc::grad(Exnerp)
          - muSponge*(U&ghat)*ghat
        )
    );

    /*volVectorField Urhs = rho*
    (
        g*thetap/thetaa
      - Cp*(thetaa + thetap)*fvc::grad(Exnerp)
    );
    volScalarField wrhs = - (Urhs & ghat);
    Urhs += wrhs*ghat;

    solve(UEqn == Urhs);
    solve(wEqn == wrhs);*/
}

