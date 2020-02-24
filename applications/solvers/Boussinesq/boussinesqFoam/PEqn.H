// Pressure equation

// Explicit part of velocity

// Rate of change of flux without pressure gradient
dVolFluxDt = 
(
    mesh.Sf() &
    (
      - fvc::interpolate(fvc::div(volFlux, u, "u"), "convection")
    )
);

// Viscous term
// Compact vector Laplacian of the velocity
Uf = linearInterpolate(u);
surfaceScalarField laplacianU
    = (fvc::interpolate(fvc::laplacian(u)) & mesh.Sf())
    + fvc::snGrad(fvc::div(volFlux - (Uf & mesh.Sf())))*mesh.magSf();

if (nu.value() > VSMALL)
{
    dVolFluxDt += nu*laplacianU;
}

// Coriolis term
if (magSqr(Omega.value()) > VSMALL)
{
    dVolFluxDt -= mesh.Sf() & (2*(Omega ^ Uf));
}

volFlux = volFlux.oldTime()
        + dt*(1-offCentre)*dVolFluxDt.oldTime()
        + dt*offCentre*dVolFluxDt;

// Pressure equation
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix PEqn
    (
        fvc::div(volFlux)
      + offCentre*dt*fvc::div(bf)
      - offCentre*dt*fvm::laplacian(P)
    );
    if (pRefCell >= 0) PEqn.setReference(pRefCell,0);
    PEqn.solve();
}

// Back-substitutions
dVolFluxDt += bf - fvc::snGrad(P)*mesh.magSf();
volFlux += offCentre*dt*(bf - fvc::snGrad(P)*mesh.magSf());
