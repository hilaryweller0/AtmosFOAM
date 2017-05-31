// Updating sigma in each partition

// Re-calculate rhoTheta (sigmaRho has been updated)
rhoTheta = sigmaRho[0]*theta[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigmaRho[ip]*theta[ip];
}

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
}

Info << "buoyant.sigma goes from " << min(sigma[1]).value() << " to "
     << max(sigma[1]).value() << endl;
