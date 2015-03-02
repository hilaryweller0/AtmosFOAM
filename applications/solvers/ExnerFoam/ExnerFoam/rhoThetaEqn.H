// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + (1-offCentre)*divPhi + offCentre*fvc::div(phi)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, theta)
  + (1-offCentre)*divPhitheta.oldTime()
  + offCentre*fvc::div(phi, theta)
);
thetaEqn.solve();

