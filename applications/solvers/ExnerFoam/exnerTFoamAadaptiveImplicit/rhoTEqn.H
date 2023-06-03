// update density according to the continuity equation
solve
(
    EulerDdt.fvmDdt(rho)
  + fvc::div(phiOld)
  + fvc::div(phiNew)
);
rhof = fvc::interpolate(rho);

// T equation
fvScalarMatrix TEqn
(
    EulerDdt.fvmDdt(rho, T)
  + fvm::div(phiNew, T, "div(phi,T)")
  + fvc::div(phiOld, T.oldTime(), "div(phi,T)")
  //- fvm::laplacian(thermophysicalTransport->alphaEff(), T)
  + thermo.p()/Cv*fvc::div(phi/rhof)
);

if (implicitT)
{
    TEqn.solve();
}
else
{
    T = TEqn.H()/TEqn.A();
}
theta = T/Exner;
thetaf = fvc::interpolate(theta);

