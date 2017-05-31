// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + (1-offCentre)*divU + offCentre*fvc::div(U)
);
