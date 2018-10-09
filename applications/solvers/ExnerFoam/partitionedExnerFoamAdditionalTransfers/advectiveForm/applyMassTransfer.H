{
// Limit and apply the mass transfer and add the mass transfer to the sum for 
// other variables

//Info << "Before limiting"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

if (alphaM == scalar(0))
{
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
}

// Transfer terms for other variables
const scalar minSigmaDiv = max(minSigma, VSMALL);
transfer[0][1] += dt*massTransfer[0][1]/max(sigmaRho[1], minSigmaDiv*rho[1]);
transfer[1][0] += dt*massTransfer[1][0]/max(sigmaRho[0], minSigmaDiv*rho[0]);
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            transfer[ip][jp] = dt*massTransfer[ip][jp]/max
            (
                sigmaRho[ip], minSigmaDiv*rho[ip]
            );
        }
    }
}

//Info << "Before mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value()
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value()
//     << "\nmin(rho[0]) = " << min(rho[0]).value() 
//     << "\nmin(rho[1]) = " << min(rho[1]).value() << endl;

// Update mass transfer terms to be either explicit or implicit depending on alphaC
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip][jp] = 1/dt * transfer[ip][jp] * 
            max
            (
                sigmaRho[ip], minSigmaDiv*rho[ip]
            ) /
            (
                1 + alphaM*(transfer[ip][jp] + transfer[jp][ip])
            );
            
        }
    }
}

for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            transfer[ip][jp] = transfer[ip][jp]*
            max
            (
                sigmaRho[ip], minSigmaDiv*rho[ip]
            ) /
            max
            (
                sigmaRho[jp], minSigmaDiv*rho[jp]
            );
        }
    }
}

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
    Info << "Min sigmaRho " << ip << ": " << min(sigmaRho[ip]).value() << endl;

    // Update the face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}



//Info << "After mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value() << endl;

}
