// update density according to the continuity equation
solve
(
    EulerDdt.fvmDdt(rho)
  + fvc::div(phiMid)
);
rhof = fvc::interpolate(rho);

