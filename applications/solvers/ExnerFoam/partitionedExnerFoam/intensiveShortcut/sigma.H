// Updating sigma in each partition

volScalarField sumSigma("sumSigma", 0*sigma[0]);

for(label ip = 0; ip < nParts; ip++)
{
    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
    sumSigma += sigma[ip];

    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];
}

Info << "1-sumSigma goes from " << min(1-sumSigma).value() << " to "
     << max(1-sumSigma).value() << endl;
Info << "buoyant.sigma goes from " << min(sigma[1]).value() << " to "
     << max(sigma[1]).value() << endl;
Info << "buoyant.rho goes from " << min(rho[1]).value() << " to "
     << max(rho[1]).value() << endl;
Info << "buoyant.theta goes from " << min(theta[1]).value() << " to "
     << max(theta[1]).value() << endl;
