// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
  + 2*rho*(Omega^(U-Ug))
  - muSponge*rho*(sqr(ghat) & U)
);

solve
(
    UEqn == rho*(g - Cp*theta*fvc::grad(Exner))
);

