// update density according to the continuity equation
drhodt = -fvc::div(U);
rho = rho.oldTime() + dt*((1-offCentre)*drhodt.oldTime() + offCentre*drhodt);


