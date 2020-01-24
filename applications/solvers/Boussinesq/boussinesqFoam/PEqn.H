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
//+ bf;

// Viscous term
// Compact vector Laplacian of the velocity
Uf = linearInterpolate(u);
//surfaceScalarField laplacianU
//    = fvc::snGrad(fvc::div(volFlux))*mesh.magSf()
//    - (linearInterpolate(fvc::curl(fvc::curl(Uf))) & mesh.Sf());

surfaceScalarField laplacianU
    = fvc::interpolate(fvc::laplacian(u)) & mesh.Sf();
laplacianU += fvc::snGrad(fvc::div(volFlux))*mesh.magSf()
            - fvc::snGrad(fvc::div(Uf & mesh.Sf()))*mesh.magSf();

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
    PEqn.setReference(0,0);
    PEqn.solve();
}

// Back-substitutions
//dVolFluxDt -= fvc::snGrad(P)*mesh.magSf();
//volFlux -= offCentre*dt*(fvc::snGrad(P)*mesh.magSf());
dVolFluxDt += bf - fvc::snGrad(P)*mesh.magSf();
volFlux += offCentre*dt*(bf - fvc::snGrad(P)*mesh.magSf());
