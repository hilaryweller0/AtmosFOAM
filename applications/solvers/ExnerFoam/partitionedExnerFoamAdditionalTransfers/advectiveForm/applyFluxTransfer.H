for(label ip = 0; ip < nParts; ip++)
{
    flux[ip] = fvc::interpolate
    (
        sigmaRhoPreTransfer[ip], "massTransfer"
    ) * volFlux[ip];
}
partitionedSurfaceScalarField fluxPreTransfer = flux;

// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            fluxTransferTerm[ip][jp] = 1/dt * transferf[ip][jp] /
            (
                1 + alphaV*(transferf[ip][jp] + transferf[jp][ip])
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
            flux[ip] += dt *
            (
                fluxPreTransfer[jp] * fluxTransferTerm[jp][ip]
              - fluxPreTransfer[ip] * fluxTransferTerm[ip][jp]
            );
        }
    }
}

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = flux[ip] / max
    (
        sigmaRhof[ip], 
        minSigmaDiv*fvc::interpolate(rho[ip], "massTransfer")
    );
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}
sigmaFlux.updateSum();
volFlux.updateSum();
