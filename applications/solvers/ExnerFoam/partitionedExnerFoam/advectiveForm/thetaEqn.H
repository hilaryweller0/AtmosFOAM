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
    volScalarField denominator
    (
        IOobject("denom", runTime.timeName(), mesh),
        mesh, dimensionedScalar("", dimless, scalar(1))
    );
    volScalarField numerator = theta[ip];

    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            volScalarField incr = dt*massTransfer[jp][ip]
                /max(sigmaRho[ip], minSigmaRho);
            numerator += incr*theta[jp];
            denominator += incr;
        }
    }
    
    theta[ip] = numerator/denominator;
}
// Update the sum
theta.updateSum();

