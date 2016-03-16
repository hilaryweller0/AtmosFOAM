// update density according to the continuity equation
rho = rho.oldTime() - dt*((1-offCentre)*divU + offCentre*fvc::div(U));


