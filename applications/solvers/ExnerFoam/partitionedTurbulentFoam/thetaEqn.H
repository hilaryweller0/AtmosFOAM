// advective form theta equation
for(label ip = 0; ip < nParts; ip++)
{
    fvScalarMatrix thetaEqn
    (
        fvm::ddt(theta[ip])
      + fvm::div(volFlux[ip], theta[ip], "theta")
      - theta[ip]*fvc::div(volFlux[ip])
      - fvm::laplacian(turbulence[ip].alphaEff()/rho[ip], theta[ip])
    );
    thetaEqn.solve();
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

theta[0] = ((1+transfer[0])*theta[0] + transfer[1]*theta[1])
           /(1+transfer.sum());
theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);

// Update the sum
theta.updateSum();

// Update the thermodynamics
for(label ip = 0; ip < nParts; ip++)
{
    thermos[ip].T() = theta[ip]*Exner;
}
