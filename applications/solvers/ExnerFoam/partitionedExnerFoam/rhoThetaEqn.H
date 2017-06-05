// flux form rho theta equation
for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] = 
    (
        theta[ip].oldTime()*sigmaRho[ip].oldTime()
      - dt*fvc::div(flux[ip], theta[ip])
    )/sigmaRho[ip];
}

// Heat in the mass transfers between partitions
for(label ip = 0; ip < nParts; ip++) { for(label jp=0; jp < nParts; jp++)
{
    if (ip != jp)
    {
        theta[ip] += dt*
        (
            massTransfer[jp][ip]*theta[jp] - massTransfer[ip][jp]*theta[ip]
        )/sigmaRho[ip];
    }
}}

rhoTheta = sigmaRho[0]*theta[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigmaRho[ip]*theta[ip];
}

// Heat transfer between partitions - a proportion of the amount to equalise
// the temperatures
volScalarField H
(
    "H",
    heatTransferCoeff*(theta[1] - theta[0])/dt*sigmaRho[0]*sigmaRho[1]/rhoSum
);

/*volScalarField Tsum = rhoTheta/rhoSum*Exner;
volScalarField H
(
    "H",
    heatTransferCoeff*thermalConductivity*(theta[1] - theta[0])/(Cp*heatDist)*pow
    (
        magg*rhoSum*Cp/(nu*Tsum*thermalConductivity)
        *mag((theta[1] - theta[0])),
        1./3
    )
);

H = min
(
    H, 
    0.25*mag(rho[1]*theta[1] - rho[0]*theta[0])/dt
);
H = max
(
    H, 
    -0.25*mag(rho[1]*theta[1] - rho[0]*theta[0])/dt
);
*/

for(label ip = 0; ip < nParts; ip++)
{
    theta[ip] += dt*H/sigmaRho[ip];
    H *= -1;
}

