{
// Limit and apply the mass transfer and add the mass transfer to the sum for 
// other variables

//Info << "Before limiting"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

// Constraint to keep sigmaRho > 0
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max
                (
                    (max(sigmaRho[ip], densityZero)-minSigma*rho[ip])/dt,
                    massTransferZero
                )
            );
        }
    }
}

// Transfer terms for other variables
const scalar minSigmaDiv = max(minSigma, VSMALL);
transfer[0] += dt*massTransfer[0][1]/max(sigmaRho[1], minSigmaDiv*rho[1]);
transfer[1] += dt*massTransfer[1][0]/max(sigmaRho[0], minSigmaDiv*rho[0]);
transfer.updateSum();

transferf[0] = fvc::interpolate(transfer[0], "massTransfer")
             + dt*dragCommon
             *fvc::interpolate(sigma[0]/rho[1], "massTransfer");
transferf[1] = fvc::interpolate(transfer[1], "massTransfer")
             + dt*dragCommon
             *fvc::interpolate(sigma[1]/rho[0], "massTransfer");
transferf.updateSum();

//Info << "Before mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value()
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value()
//     << "\nmin(rho[0]) = " << min(rho[0]).value() 
//     << "\nmin(rho[1]) = " << min(rho[1]).value() << endl;

// Apply mass transfer terms (operator split) to sigmaRho equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigmaRho[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
        }
    }

    // Update the face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}

//Info << "After mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value() << endl;

theta[0] = ((1+transfer[0])*theta[0] + transfer[1]*theta[1])
           /(1+transfer.sum());
theta[1] = (theta[1] + transfer[0]*theta[0])/(1+transfer[0]);

thetaVar[0] = ((1+transfer[0])*thetaVar[0] + transfer[1]*thetaVar[1])
           /(1+transfer.sum());
thetaVar[1] = (thetaVar[1] + transfer[0]*thetaVar[0])/(1+transfer[0]);

volFlux[0] =
(
    (1+transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
)/(1+transferf.sum());
volFlux[1] = (volFlux[1] + transferf[0]*volFlux[0])
            /(1+transferf[0]);


// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}

sigmaRho.updateSum();
theta.updateSum();
thetaVar.updateSum();
sigmaFlux.updateSum();
volFlux.updateSum();

}
