// update density according to the continuity equation
solve(fvm::ddt(rho) + fvc::div(U));


