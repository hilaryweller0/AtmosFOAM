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
    EulerDdt.fvmDdt(rho, thetap)
  + fvm::div(phiNew, thetap, "div(phi,thetap)")
  + fvc::div(phiOld, thetap.oldTime(), "div(phi,thetap)")
  + (1-ocAlpha)*rho.oldTime()*w.oldTime()*dThetaadz
  + ocAlpha*rho*w*dThetaadz
  //- fvm::laplacian(thermophysicalTransport->alphaEff(), theta)
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
thetaf = thetapf + thetaaf;
theta = thetaa + thetap;
