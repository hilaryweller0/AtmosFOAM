// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] = 
    (
        theta[ip].oldTime()*sigma[ip]*rho[ip].oldTime()
      - dt*fvc::div(flux[ip], sigma[ip]*theta[ip])
    )/(sigma[ip]*rho[ip]);
}

// Heat in the mass transfers between partitions
for(label ip = 0; ip < nParts; ip++) { for(label jp=0; jp < nParts; jp++)
{
    if (ip != jp)
    {
        theta[ip] += dt*
        (
            massTransfer[jp][ip]*theta[jp] - massTransfer[ip][jp]*theta[ip]
        )/(sigma[ip]*rho[ip]);
    }
}}

rhoTheta = sigma[0]*rho[0]*theta[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigma[ip]*rho[ip]*theta[ip];
}

// Heat transfer between partitions - a proportion of the amount to equalise
// the temperatures
volScalarField H
(
    "H",
    heatTransferCoeff*(theta[1] - theta[0])/dt*sigma[0]*rho[0]*sigma[1]*rho[1]
    /rhoSum
);

for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] += dt*H/(sigma[ip]*rho[ip]);
    H *= -1;
}

