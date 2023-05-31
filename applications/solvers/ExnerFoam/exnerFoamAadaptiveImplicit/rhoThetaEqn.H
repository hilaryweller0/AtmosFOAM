// update density according to the continuity equation
solve
(
    EulerDdt.fvmDdt(rho)
  + fvc::div(phiOld)
  + fvc::div(phiNew)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, theta)
  + fvm::div(phiNew, theta, "div(phi,theta)")
  + fvc::div(phiOld, theta.oldTime(), "div(phi,theta)")
  //- fvm::laplacian(thermophysicalTransport->alphaEff(), theta)
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
