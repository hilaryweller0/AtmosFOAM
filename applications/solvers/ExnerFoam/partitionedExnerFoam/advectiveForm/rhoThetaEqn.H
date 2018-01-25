// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] =-fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);

    // Include the mass transfer terms
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            theta.ddt()[ip] += massTransfer[jp][ip]/sigmaRho[ip]*
                               (theta[jp] - theta[ip]);
        }
    }
    
    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*
    (
        (1-offCentre)*theta.ddt()[ip].oldTime()
      + offCentre*theta.ddt()[ip]
    );
}
// Update the sum
theta.updateSum();

