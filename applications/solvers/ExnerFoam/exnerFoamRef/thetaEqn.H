// theta equation
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + fvc::div(phiMid, thetap, "div(phi,thetap)")
  + rho*w*(fvc::grad(offCentre*thetaaf) & kdir)
  - fvm::laplacian
    (
        offCentre*fvc::interpolate(thermophysicalTransport->alphaEff()),
        thetap
    )
  - thetapRHS
);

thetaEqn.solve("thetaEqn");

thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
