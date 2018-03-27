// advective form theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] =-fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);

    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*
    (
        (1-offCentre)*theta.ddt()[ip].oldTime()
      + offCentre*theta.ddt()[ip]
    );
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

//// Explicit application of mass transfer terms using conserved variables
//theta[0] += dt*(massTransfer[1][0]*theta[1] - massTransfer[0][1]*theta[0])
//            / max(sigmaRho[0], minSigmaDiv*rho[0]);
//theta[1] += dt*(massTransfer[0][1]*theta[0] - massTransfer[1][0]*theta[1])
//            / max(sigmaRho[1], minSigmaDiv*rho[1]);

// Update the sum
theta.updateSum();

