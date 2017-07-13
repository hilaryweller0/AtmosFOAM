// Momentum predictor
fvVectorMatrix UEqn
(
    fvm::ddt(rho,U)
  + fvc::div(phi, U)
  + (2*rho*(Omega^(U-Ug)))
);

if (piso.momentumPredictor())
{
    solve(UEqn == - Cp*rho*theta*fvc::grad(Exner) + rho*g);
}
