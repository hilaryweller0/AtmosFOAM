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
  + offCentre*fvm::div(phi, theta)
  - fvm::laplacian(turbulence->alphaEff(), theta)
);
// note: alphaEff() would take into account the effects of molecular viscosity

thetaEqn.solve();

thermo.T() = theta*Exner;
//thermo.he() == thermo.he(thermo.p(),thermo.T());
//thermo.correct();
