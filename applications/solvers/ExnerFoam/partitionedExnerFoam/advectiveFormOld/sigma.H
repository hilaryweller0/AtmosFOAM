// Updating sigma in each partition

volScalarField sumSigma("sumSigma", 0*sigma[0]);

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
    sumSigma += sigma[ip];
    
    Info << sigma[ip].name() << " goes from " << min(sigma[ip]).value()
         << " to " << max(sigma[ip]).value() << endl;
}
Info << "1-sumSigma goes from " << min(1-sumSigma).value()
     << " to " << max(1-sumSigma).value() << endl;

