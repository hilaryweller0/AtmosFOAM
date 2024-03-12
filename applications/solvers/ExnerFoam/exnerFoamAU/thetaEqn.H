// theta equation
alphaEff = fvc::interpolate(thermophysicalTransport->kappaEff());
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, theta)
  + fvm::div(phi, theta)
  - fvm::laplacian(alphaEff, theta)
);

if (implicitT)
{
    thetaEqn.solve();
}
else
{
    theta = thetaEqn.H()/thetaEqn.A();
}
thetaf = fvc::interpolate(theta);
