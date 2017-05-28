// update rhoSigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRho[ip] = sigmaRho[ip].oldTime() - dt*fvc::div(flux[ip]);
}

