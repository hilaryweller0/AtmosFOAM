// update density according to the continuity equation

dRhodt = -fvc::div(phi);

solve
(
    EulerDdt.fvmDdt(rho)
  - ocAlpha*dRhodt
  - (1-ocAlpha)*dRhodt.oldTime()
);
rhof = fvc::interpolate(rho);

