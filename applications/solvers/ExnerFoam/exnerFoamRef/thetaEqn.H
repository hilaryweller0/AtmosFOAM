// theta equation
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + fvm::div(aphi, thetap, "div(phi,thetap)")
  + rho*w*(fvc::grad(offCentre*thetaaf) & kdir)
  - fvm::laplacian
    (
        offCentre*fvc::interpolate(thermophysicalTransport->alphaEff()),
        thetap
    )
  - thetapRHS
);

if (implicitT)
{
    thetaEqn.solve();
}
else
{
    thetap = thetaEqn.H()/thetaEqn.A();
}

thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
