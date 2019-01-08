// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(volFlux, U)
  - fvm::Sp(fvc::div(volFlux), U)
  + turbulence->divDevRhoReff(U)
  + 2*(Omega^(U-Ug))
  - muSponge*(sqr(ghat) & U)
);

solve
(
    UEqn ==
    fvc::reconstruct(gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf())
);

Uf = fvc::interpolate(U);

