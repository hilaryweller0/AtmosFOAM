// Updating sigma in each partition

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = theta.sum()/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
    
    Info << sigma[ip].name() << " goes from " << min(sigma[ip]).value()
         << " to " << max(sigma[ip]).value() << endl;
}

rho.updateSum();
sigma.updateSum();

Info << "1-sigma.sum() goes from " << min(1-sigma.sum()).value()
     << " to " << max(1-sigma.sum()).value() << endl;

