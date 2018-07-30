// Momentum predictor
fvVectorMatrix UEqn
(
    fvm::ddt(rho,U)
  + fvc::div(phi, U)
  + (2*rho*(Omega^(U-Ug)))
  - nu*fvm::laplacian(rho, U)
  + fvm::Sp(rho*muSponge, U)
  - rho*muSponge*Urelax
);

if (pimple.momentumPredictor())
{
    solve(UEqn == rho*(g - Cp*theta*fvc::grad(Exner)));
}
