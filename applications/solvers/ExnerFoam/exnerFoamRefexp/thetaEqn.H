// theta equation

dRhoThetadt = -fvc::div(phi, thetap)
              - rho*U.component(2)*dThetaadz;

fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  - ocAlpha*dRhoThetadt
  - (1-ocAlpha)*dRhoThetadt.oldTime()
);

thetaEqn.solve();

thetapf = fvc::interpolate(thetap);
thetaf = thetapf + thetaaf;
theta = thetaa + thetap;
