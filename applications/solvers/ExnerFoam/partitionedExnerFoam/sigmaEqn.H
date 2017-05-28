{
// update sigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] = (sigmaRho[ip].oldTime() - dt*fvc::div(flux[ip]))/rho[ip];
}

// Check that sum(sigma) = 1
volScalarField sumSigma = sigma[0];
for(label ip = 1; ip < nParts; ip++)
{
    sumSigma += sigma[ip];
}
Info << "1-sumSigma goes from " << 1-max(sumSigma).value() << " to "
     << 1-min(sumSigma).value() << endl;
}
