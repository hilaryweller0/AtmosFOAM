// Updating rho in each partition

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];
}

Info << "buoyant.rho goes from " << min(rho[1]).value() << " to "
     << max(rho[1]).value() << endl;

