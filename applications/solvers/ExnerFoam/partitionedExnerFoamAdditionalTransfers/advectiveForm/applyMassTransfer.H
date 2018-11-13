{
// Limit and apply the mass transfer and add the mass transfer to the sum for 
// other variables

partitionedVolScalarField sigmaRhoPreTransfer = sigmaRho;
partitionedSurfaceScalarField sigmaRhofPreTransfer = sigmaRhof;

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
            transferf[ip][jp] = fvc::interpolate(transfer[ip][jp], "massTransfer");
        }
    }
}

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



Info << "After mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value() << endl;

//Conditions for the treatment of sigmaRho in the transfer terms.
partitionedVolScalarField sigmaRhoNumerator = sigmaRho;
partitionedVolScalarField sigmaRhoDenominator = sigmaRho;
partitionedSurfaceScalarField sigmaRhofNumerator = sigmaRhof;
partitionedSurfaceScalarField sigmaRhofDenominator = sigmaRhof;
if (transferNumeratorOld)
{
    sigmaRhoNumerator = sigmaRhoPreTransfer;
}
if (transferDenominatorOld)
{
    sigmaRhoDenominator = sigmaRhoPreTransfer;
}

//Ensure denominator not zero
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRhoNumerator[ip] = max(sigmaRhoNumerator[ip], minSigmaDiv*rho[ip]);
    sigmaRhoDenominator[ip] = max(sigmaRhoDenominator[ip], minSigmaDiv*rho[ip]);
    sigmaRhofNumerator[ip] = fvc::interpolate
    (
        sigmaRhoNumerator[ip], "massTransfer"
    );
    sigmaRhofDenominator[ip] = fvc::interpolate
    (
        sigmaRhoDenominator[ip], "massTransfer"
    );
}

//Transfer theta 
if (useAdvectiveThetaTransfer)
{
    #include "applyThetaTransfer.H"
}
else if (useFluxThetaTransfer)
{
    #include "applyThetaFluxTransfer.H"
}

//Transfer velocity
if (useVelocityTransfer)
{
    #include "applyVelocityTransfer.H"
}
else if (useVolFluxTransfer)
{
    #include "applyVolFluxTransfer.H"
}
else if (useMomentumTransfer)
{
    #include "applyMomentumTransfer.H"
}
else if (useFluxTransfer)
{
    #include "applyFluxTransfer.H"
}

}
