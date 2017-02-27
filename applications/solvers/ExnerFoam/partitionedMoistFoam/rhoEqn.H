// Advect the total density
dRhodt = -fvc::div(flux);
rho = rho.oldTime() + dt*
(
    (1-offCentre)*dRhodt.oldTime() + offCentre*dRhodt
);


