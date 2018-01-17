// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] = theta[ip].oldTime() - dt*
    (
        fvc::div(volFlux[ip], theta[ip]) - theta[ip]*fvc::div(volFlux[ip])
    );
}

rhoTheta = sigmaRho[0]*theta[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigmaRho[ip]*theta[ip];
}

