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
  //- fvm::laplacian(thermophysicalTransport->alphaEff(), theta)
);

if (implicitT)
{
    thetaEqn.solve();
}
theta = thetaEqn.H()/thetaEqn.A();
thetaf = fvc::interpolate(theta);

surfaceScalarField gradTheta("gradTheta", fvc::snGrad(theta));
surfaceScalarField ghGradThetaByTheta
(
    "ghGradThetaByTheta",
    ghf*gradTheta/thetaf
);
