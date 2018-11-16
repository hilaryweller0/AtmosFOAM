for(label ip = 0; ip < nParts; ip++)
{
    massFlux[ip] = fvc::interpolate
    (
        sigmaRhoPreTransfer[ip], "massTransfer"
    ) * volFlux[ip];
}
partitionedSurfaceScalarField massFluxPreTransfer = massFlux;

// Set flux transfer terms to be either explicit or implicit based on alphaV
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massFluxTransferTerm[ip][jp] = 1/dt * transferf[ip][jp] /
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
            massFlux[ip] += dt *
            (
                massFluxPreTransfer[jp] * massFluxTransferTerm[jp][ip]
              - massFluxPreTransfer[ip] * massFluxTransferTerm[ip][jp]
            );
        }
    }
}

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    volFlux[ip] = massFlux[ip] / max
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
