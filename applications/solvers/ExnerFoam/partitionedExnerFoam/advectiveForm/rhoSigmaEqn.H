// update rhoSigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRho[ip] = sigmaRho[ip].oldTime()
                 - dt*fvc::div(volFlux[ip], sigmaRho[ip]);
}

rhoSum = sigmaRho[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoSum += sigmaRho[ip];
}

Info << "buoyant sigma goes from " << min(sigmaRho[1]/rho[1]).value()
     << " to " << max(sigmaRho[1]/rho[1]).value() << endl;

