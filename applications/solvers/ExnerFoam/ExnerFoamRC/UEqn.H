// Momentum predictor
fvVectorMatrix UEqn
(
    fvm::ddt(rho,U)
  + fvc::div(phi, U)
  + (2*rho*(Omega^(U-Ug)))
  - nu*fvm::laplacian(rho, U)
  + rho/runTime.deltaT()*(muSponge&U)
);

if (pimple.momentumPredictor())
{
    solve(UEqn == rho*(g - Cp*theta*fvc::grad(Exner)));
}
