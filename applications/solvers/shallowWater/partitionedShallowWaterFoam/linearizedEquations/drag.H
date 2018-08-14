//Calculate drag.
drag[ip] *= 0;
if (useDrag)
{
    if (0)
    {
        //volScalarField numberOfBubbles = 4*(Nmax-1)*sigma[ip]*(1-sigma[ip])*Foam::exp(-Nmax*lengthScale*mag(fvc::grad(sigma[ip]))) + scalar(1);
        volScalarField numberOfBubbles = 4*(Nmax-1)*sigma[ip]*(1-sigma[ip])*(scalar(1)-lengthScale*mag(fvc::grad(sigma[ip]))) + scalar(1);
        //volScalarField bubbleRadius = max(bubbleRadiusMin,0.5*lengthScale*sigma[ip]/numberOfBubbles);
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                volScalarField bubbleRadius = max(bubbleRadiusMin,bubbleRadiusMax*sigma[ip]*sigma[jp]);
                drag[ip] -= 0.5*fvc::interpolate(sigma[jp])*dragCoefficient/fvc::interpolate(bubbleRadius) * mag(Uf[ip]-Uf[jp]) * (Uf[ip]-Uf[jp]);
            }
        }
    }
    
    if (1)
    {
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                volScalarField bubbleRadius = max(bubbleRadiusMin,bubbleRadiusMax*sigma[ip]*sigma[jp]);
                surfaceScalarField bubbleRadiusf = fvc::interpolate(bubbleRadius);
                surfaceScalarField sigmajf = fvc::interpolate(sigma[jp]);
                forAll(drag[ip],celli)
                {
                    if ( magSqr(0.5*dragCoefficient/bubbleRadiusf[celli] * mag(Uf[ip][celli]-Uf[jp][celli])) < 1)
                    {
                        drag[ip][celli] -= 0.5*sigmajf[celli]*dragCoefficient/bubbleRadiusf[celli] * mag(Uf[ip][celli]-Uf[jp][celli]) * (Uf[ip][celli]-Uf[jp][celli]);
                    }
                    else
                    {
                        drag[ip][celli] -= sigmajf[celli] * (Uf[ip][celli]-Uf[jp][celli]);
                    }
                }
            }
        }
    }
}
