partitionedSurfaceScalarField volFluxPreTransfer = volFlux;

// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            volFluxTransferTerm[ip][jp] = 1/dt * transferf[ip][jp] *
            (
                sigmaRhofNumerator[ip] / sigmaRhofDenominator[jp]
            ) /
            (
                1 + alphaV *
                (
                    transferf[ip][jp]  *
                    (
                        sigmaRhofNumerator[ip] / sigmaRhofDenominator[jp]
                    )
                  + transferf[jp][ip]  *
                    (
                        sigmaRhofNumerator[jp] / sigmaRhofDenominator[ip]
                    )
                )
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
            volFlux[ip] += dt * volFluxTransferTerm[jp][ip] * 
            (
                volFluxPreTransfer[jp]
              - volFluxPreTransfer[ip]
            );
        }
    }
}

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}
sigmaFlux.updateSum();
volFlux.updateSum();
