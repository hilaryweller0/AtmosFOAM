
rhof = fvc::interpolate(rho);
thetaRho = fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl), "theta");

// Rate of change of flux
dVdt = 
     -(H.delta() &
     (
        fvc::interpolate(fvc::div(U,u))
      + (2*rhof*(Omega^(Uf-Ug)))
     ))
    + rhof*(gd - Cp*thetaRho*fvc::snGrad(Exner)*H.magd());

// Update V (momenum in the d direction) and U (normal mass flux)
V = V.oldTime() + offCentre*dt*dVdt;
U = H.ddirToFlux(V);
u = H.reconstructd(V)/rho;
Uf = fvc::interpolate(u);

