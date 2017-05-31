// update rhoSigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRho[ip] = sigmaRho[ip].oldTime() - dt*fvc::div(flux[ip]);
}

rhoSum = sigmaRho[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoSum += sigmaRho[ip];
}

