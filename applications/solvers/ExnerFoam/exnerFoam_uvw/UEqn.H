// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
);

fvScalarMatrix wEqn(fvm::Sp(rho*muSponge, w));

if (SIgravityWaves)
{
    volScalarField G = -ocAlpha*runTime.deltaT()*Cp*rho
         *(fvc::grad(Exner) & (-ghat)) * (fvc::grad(theta) & (-ghat));
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));

    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix(fvm::ddt(rho, w) + fvm::div(phi, w));
}

if (impU)
{
    volVectorField Urhs = rho*fvc::reconstruct
    (
        gSf
      - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
    );
    volScalarField wrhs = - (Urhs & ghat);
    Urhs += wrhs*ghat;

    solve(UEqn == Urhs);
    solve(wEqn == wrhs);
}
