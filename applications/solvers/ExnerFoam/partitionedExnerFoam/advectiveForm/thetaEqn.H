// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta.ddt()[ip] =-fvc::div(volFlux[ip], theta[ip], "theta")
                    + theta[ip]*fvc::div(volFlux[ip]);

    // Update theta according to the rate of change
    theta[ip] = theta[ip].oldTime() + dt*
    (
        (1-offCentre)*theta.ddt()[ip].oldTime()
      + offCentre*theta.ddt()[ip]
    );
}

// Mass transfer terms
for(label ip = 0; ip < nParts; ip++)
{
// Include the mass transfer terms (operator split, explicit source)
//    for(label jp=0; jp < nParts; jp++)
//    {
//        if (ip != jp)
//        {
//            theta[ip] += dt*massTransfer[jp][ip]/sigmaRho[ip]*
//                         (theta[jp] - theta[ip]);
//        }
//    }

    // Include the mass transfer terms (operator split, implicit source)
    volScalarField denominator = sigmaRho[ip];
    volScalarField numerator = sigmaRho[ip]*theta[ip];

    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            numerator += dt*massTransfer[jp][ip]*theta[jp];
            denominator += dt*massTransfer[jp][ip];
        }
    }
    
    for(label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (mag(denominator[cellI]) > 1e-4)
        {
            theta[ip][cellI] = numerator[cellI]/denominator[cellI];
        }
    }
}
// Update the sum
theta.updateSum();

