// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

solve
(
    UEqn ==
    fvc::reconstruct(rhof*(gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()))
);

