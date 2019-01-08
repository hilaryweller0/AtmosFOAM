// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + fvc::div(phi)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(theta)
  + fvm::div(volFlux, theta)
  - fvm::Sp(fvc::div(volFlux), theta)
  - fvm::laplacian(turbulence->alphaEff()/rho, theta)
);

theta.correctBoundaryConditions();
thetaEqn.solve();
thetaf = fvc::interpolate(theta);

thermo.T() = theta*Exner;
//thermo.he() == thermo.he(thermo.p(),thermo.T());
//thermo.correct();
