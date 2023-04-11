// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + fvc::div(phi)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, theta)
  + fvm::div(phi, theta)
  - fvm::laplacian(thermophysicalTransport->alphaEff(), theta)
);

thetaEqn.solve();
/*// Explitic solve for energy conservation
solve
(
    fvm::ddt(rho, theta)
  + fvc::div(phi, theta)
  - fvc::laplacian(thermophysicalTransport->alphaEff(), theta)
);
*/
