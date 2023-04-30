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
);

if (implicitT) thetaEqn.solve();
else theta = thetaEqn.H()/thetaEqn.A();

thetaf = fvc::interpolate(theta);
