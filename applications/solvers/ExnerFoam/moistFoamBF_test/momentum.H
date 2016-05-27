thetaRho = fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl), "theta");

// Rate of change of flux
dVdt = -(H.delta() & fvc::interpolate(fvc::div(U,u) - u*divu))
      + gd
      - Cp*thetaRho*fvc::snGrad(Exner)*H.magd();

// Update V (momenum in the d direction) and U (normal mass flux)
V = V.oldTime() + (1-offCentre)*dt*dVdt.oldTime() + offCentre*dt*dVdt;
U = H.ddirToFlux(V);
u = H.reconstructd(V);
Uf = fvc::interpolate(u);

