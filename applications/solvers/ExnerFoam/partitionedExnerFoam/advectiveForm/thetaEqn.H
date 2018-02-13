// flux form rho theta equation
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

//// The heat transfer terms
//volScalarField heatTransfer = -thetaTransferDiffusivity*
//(
//    min
//    (
//        fvc::laplacian(theta[0], "theta"), heatTransferZero
//    )
//  - max
//    (
//        fvc::laplacian(theta[1], "theta"), heatTransferZero
//    )
//);
//theta[0] -= dt*heatTransfer;
//theta[1] += dt*heatTransfer;

// Update the sum
theta.updateSum();

