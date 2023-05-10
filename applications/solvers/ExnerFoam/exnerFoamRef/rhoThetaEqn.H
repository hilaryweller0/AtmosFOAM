// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + fvc::div(phi)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, thetap)
  + fvm::div(phi, thetap)
  - fvm::laplacian(thermophysicalTransport->alphaEff(), thetap)
  - rho*(U & ghat)*dThetadz
);

if (implicitT)
{
    thetaEqn.solve();
}
thetap = thetaEqn.H()/thetaEqn.A();
thetaf = fvc::interpolate(thetaa + thetap);
