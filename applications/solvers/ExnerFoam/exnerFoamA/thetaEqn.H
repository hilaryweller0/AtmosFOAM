// theta equation
alphaEff = fvc::interpolate(thermophysicalTransport->kappaEff());
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, theta)
  + fvm::div(ImEx*aphi, theta, "div(phi,theta)")
  + fvc::div((1-ImEx)*aphi, theta, "div(phi,theta)")
  - fvm::laplacian
    (
        offCentre*alphaEff,
        theta
    )
  - thetaRHS
);

thetaEqn.solve();
thetaf = fvc::interpolate(theta);
