// update rhoSigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Calculate the rate of change
    sigmaRho.ddt()[ip] = -fvc::div(volFlux[ip], sigmaRho[ip], "sigmaRho");

    // Update sigmaRho based on the rate of change
    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*
    (
        (1-offCentre)*sigmaRho.ddt()[ip].oldTime()
      + offCentre*sigmaRho.ddt()[ip]
    );
}

// Calculate mass tranfer terms as diffusion between partitions
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
       massTransfer[ip][jp] = 0.25*Ksigma*
        (
            rho[jp]*fvc::laplacian(sigma[jp])
          - rho[ip]*fvc::laplacian(sigma[ip])
        );
    massTransfer[jp][ip] = -massTransfer[ip][jp];
    }
}

// Stability constraints for mass transfer terms
for(label it = 0; it < 4; it++) for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Constraint to keep sigmaRho > 0
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                massTransfer[jp][ip] + sigmaRho[ip]/dt - minSigmaRho/dt
            );

            // Constraint for transported variables
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max(sigmaRho[jp]/dt, minSigmaRho/dt)
            );
        }
    }
}

Info << "Before mass transfers, sigmaRho[0] goes from "
     << min(sigmaRho[0]).value() << " to " << max(sigmaRho[0]).value() << endl;
Info << "Before mass transfers, sigmaRho[1] goes from "
     << min(sigmaRho[1]).value() << " to " << max(sigmaRho[1]).value() << endl;

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
}

Info << "After mass transfers, sigmaRho[0] goes from "
     << min(sigmaRho[0]).value() << " to " << max(sigmaRho[0]).value() << endl;
Info << "After mass transfers, sigmaRho[1] goes from "
     << min(sigmaRho[1]).value() << " to " << max(sigmaRho[1]).value() << endl;

// Update the sums and face values
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRho.updateSum();
sigmaRhof.updateSum();


