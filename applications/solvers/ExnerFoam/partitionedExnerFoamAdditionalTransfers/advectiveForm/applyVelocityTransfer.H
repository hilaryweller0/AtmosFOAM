partitionedVolVectorField velocityPreTransfer = u;

// Set theta transfer terms to be either explicit or implicit based on alphaT
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            velocityTransferTerm[ip][jp] = 1/dt * transfer[ip][jp] *
            (
                sigmaRhoNumerator[ip] / sigmaRhoDenominator[jp]
            ) /
            (
                1 + alphaV *
                (
                    transfer[ip][jp] * 
                    (
                        sigmaRhoNumerator[ip] / sigmaRhoDenominator[jp]
                    )
                  + transfer[jp][ip] * 
                    (
                        sigmaRhoNumerator[jp] / sigmaRhoDenominator[ip]
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
            u[ip] += dt * velocityTransferTerm[jp][ip] * 
            (
                velocityPreTransfer[jp]
              - velocityPreTransfer[ip]
            );
        }
    }
}

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
    volFlux[ip] = (Uf[ip] & mesh.Sf());
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
}
sigmaFlux.updateSum();
volFlux.updateSum();
