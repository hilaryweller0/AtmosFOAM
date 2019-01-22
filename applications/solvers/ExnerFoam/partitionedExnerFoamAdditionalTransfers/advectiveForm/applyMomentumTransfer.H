for(label ip = 0; ip < nParts; ip++)
{
    momentum[ip] = sigmaRhoPreTransfer[ip]*u[ip];
}
partitionedVolVectorField momentumPreTransfer = momentum;

// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            momentumTransferTerm[ip][jp] = 1/dt * transfer[ip][jp] /
            (
                1 + alphaV*(transfer[ip][jp] + transfer[jp][ip])
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
            momentum[ip] += dt *  
            (
                momentumPreTransfer[jp] * momentumTransferTerm[jp][ip]
              - momentumPreTransfer[ip] * momentumTransferTerm[ip][jp]
            );
        }
    }
}

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    u[ip] = momentum[ip] / max
    (
        sigmaRho[ip], 
        minSigmaDiv*rho[ip]
    );
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
    volFlux[ip] = (Uf[ip] & mesh.Sf());
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
}
sigmaFlux.updateSum();
volFlux.updateSum();
