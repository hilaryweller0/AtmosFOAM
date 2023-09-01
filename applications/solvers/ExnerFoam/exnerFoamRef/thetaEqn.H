// theta equation

fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + fvm::div(ImEx*aphi, thetap, "div(phi,thetap)")
  + fvc::div((1-ImEx)*aphi, thetap, "div(phi,thetap)")
  + fvc::average(offCentrePG)*rho*w*dThetaadz
  - fvm::laplacian
    (
        offCentre*fvc::interpolate(thermophysicalTransport->alphaEff()),
        thetap
    )
  - thetapRHS
);

thetaEqn.solve();

thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
