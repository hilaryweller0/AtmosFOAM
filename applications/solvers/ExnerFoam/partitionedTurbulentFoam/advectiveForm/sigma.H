// Updating sigma in each partition and sigma on the face

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = theta.sum()/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];


    // Update sigma on the face and limit to between zero and one
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
    sigmaf[ip] = min(max(sigmaf[ip], scalar(0)), scalar(1));
}

rho.updateSum();
sigma.updateSum();
sigmaf.updateSum();

