// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + fvc::div(phi)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, thetap)
  + fvm::div(phi, thetap)
  + rho*w*(fvc::grad(thetaa) & kdir)//dThetaadz
);

if (implicitT)
{
    thetaEqn.solve();
}
thetap = thetaEqn.H()/thetaEqn.A();
thetapf = fvc::interpolate(thetap);
theta = thetaa + thetap;
thetaf = thetapf + thetaaf;

