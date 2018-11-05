partitionedSurfaceScalarField volFluxPreTransfer = volFlux;

// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            surfaceScalarField transferfIJ = fvc::interpolate
            (
                transfer[ip][jp] *
                (
                    sigmaRhoNumerator[ip] / sigmaRhoDenominator[jp]
                ),
                "massTransfer"
            );
            surfaceScalarField transferfJI = fvc::interpolate
            (
                transfer[jp][ip] *
                (
                    sigmaRhoNumerator[jp] / sigmaRhoDenominator[ip]
                ),
                "massTransfer"
            );
        
            volFluxTransferTerm[ip][jp] = 1/dt * transferfIJ /
            (
                1 + alphaV *
                (
                    transferfIJ
                  + transferfJI
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
