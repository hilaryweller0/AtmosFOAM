// theta equation
alphaEff = fvc::interpolate(thermophysicalTransport->kappaEff());
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + fvm::div(ImEx*aphi, thetap, "div(phi,thetap)")
  + fvc::div((1-ImEx)*aphi, thetap, "div(phi,thetap)")
  + fvc::average(offCentreGW)*rho*w*dThetaadz
  - fvm::laplacian
    (
        offCentre*alphaEff,
        thetap
    )
  - thetapRHS
);

thetaEqn.solve();

thetapf = fvc::interpolate(thetap);
thetaf = thetapf + thetaaf;
theta = thetaa + thetap;
