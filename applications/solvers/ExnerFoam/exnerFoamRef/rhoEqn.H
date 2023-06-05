// update density according to the continuity equation
solve
(
    EulerDdt.fvmDdt(rho)
  - rhoRHS.oldTime()
  + fvc::div(aphi)
);
rhof = fvc::interpolate(rho);

