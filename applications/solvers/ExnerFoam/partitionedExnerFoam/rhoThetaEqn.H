// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] = 
    (
        theta[ip].oldTime()*sigmaRho[ip].oldTime()
      - dt*fvc::div(flux[ip], theta[ip])
    )/sigmaRho[ip];
}

