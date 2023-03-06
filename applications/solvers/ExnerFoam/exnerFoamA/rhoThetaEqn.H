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
  - fvm::laplacian(turbulence->alphaEff(), theta)
);

thetaEqn.solve();

solve
(
    fvm::ddt(rho, theta)
  + fvc::div(phi, theta)
  - fvc::laplacian(turbulence->alphaEff(), theta)
);

thetaf = fvc::interpolate(theta);

