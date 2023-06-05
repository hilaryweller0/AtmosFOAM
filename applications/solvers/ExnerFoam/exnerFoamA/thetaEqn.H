// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, theta)
  + fvm::div(phi, theta)
  - fvm::laplacian(thermophysicalTransport->alphaEff(), theta)
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
