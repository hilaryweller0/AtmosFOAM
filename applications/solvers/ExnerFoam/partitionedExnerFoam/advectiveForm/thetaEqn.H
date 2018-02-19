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

// Heat transfers
if (nParts != 2)
{
    FatalErrorIn("thetaEqn.H")
        << " heat transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}
// Mass weighted heat transfer
volScalarField heatTransfer
(
    "heatTransfer",
    -sigmaRho[0]*thetaTransferDiffusivity*min
    (
        fvc::laplacian(theta[0]),
        dimensionedScalar("", dimensionSet(0,-2,0,1,0), scalar(0))
    )
    -sigmaRho[1]*thetaTransferDiffusivity*max
    (
        fvc::laplacian(theta[1]),
        dimensionedScalar("", dimensionSet(0,-2,0,1,0), scalar(0))
    )
);
// Limit the heat transfer
heatTransfer = min
(
    heatTransfer,
    1e-3*min(sigmaRho[0]*theta[0], sigmaRho[1]*theta[1])/dt
);
heatTransfer = max
(
    heatTransfer,
   -1e-3*min(sigmaRho[0]*theta[0], sigmaRho[1]*theta[1])/dt
);
heatTransfer.write();

theta[0] -= dt*heatTransfer/max(sigmaRho[0], minSigma*rho[0]);
theta[1] += dt*heatTransfer/max(sigmaRho[1], minSigma*rho[1]);

// Update the sum
theta.updateSum();

