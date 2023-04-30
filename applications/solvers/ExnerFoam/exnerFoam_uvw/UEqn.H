// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
);

// Implicit solution of w equation
fvScalarMatrix wEqn(fvm::Sp(rho*muSponge, w));

if (SIgravityWaves || hydrostatic)
{
    surfaceScalarField Gf = -sqr(rhof)*Cp*fvc::snGrad(theta)*fvc::snGrad(Exner)
                            /fvc::interpolate(thetaEqn.A());

    // Apply this only where on vertical pointing faces where G>0
    surfaceScalarField applySIg = mag(gUnitNormal)*max(sign(Gf), scalar(0));
    
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());

    volScalarField G = fvc::average(Gf);
    G = max(G, dimensionedScalar("", G.dimensions(), scalar(0)));
    
    wEqn += fvScalarMatrix(fvm::Sp(G, w));
}

if (!hydrostatic)
{
    wEqn += fvScalarMatrix(fvm::ddt(rho, w) + fvm::div(phi, w));
}

if (implicitU)
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
