// update h and sigma in each partition and calculate hSum

// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    hf[ip] = fvc::interpolate(h[ip]);
    h[ip] = h[ip].oldTime() - dt*fvc::div(volFlux[ip]*hf[ip]);

    if (ip == 0)
    {
        hSum = h[ip];
        hfSum = hf[ip];
    }
    else
    {
        hSum += h[ip];
        hfSum += hf[ip];
    }
}
// Update sigma and hScale
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] = h[ip]/hSum;
    Info << "sigma[" << ip << "] goes from "
         << min(sigma[ip]).value() << " to "
         << max(sigma[ip]).value() << endl;
}

