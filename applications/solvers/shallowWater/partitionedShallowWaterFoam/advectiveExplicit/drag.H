//Calculate drag.
if (modelDrag)
{
    if (1)
    {
        drag[ip] *= 0;
        //volScalarField numberOfBubbles = 4*(Nmax-1)*sigma[ip]*(1-sigma[ip])*Foam::exp(-Nmax*lengthScale*mag(fvc::grad(sigma[ip]))) + scalar(1);
        volScalarField numberOfBubbles = 4*(Nmax-1)*sigma[ip]*(1-sigma[ip])*(scalar(1)-lengthScale*mag(fvc::grad(sigma[ip]))) + scalar(1);
        //volScalarField bubbleRadius = max(bubbleRadiusMin,0.5*lengthScale*sigma[ip]/numberOfBubbles);
        for (label jp = 0; jp < ip; jp++)
        {
            volScalarField bubbleRadius = max(bubbleRadiusMin,bubbleRadiusMax*sigma[ip]*sigma[jp]);
            drag[ip] -= 0.5*fvc::interpolate(sigma[jp])*dragCoefficient/fvc::interpolate(bubbleRadius) * mag(Uf[ip]-Uf[jp]) * (Uf[ip]-Uf[jp]);
        }
        for (label jp = ip+1; jp < nParts; jp++)
        {
            volScalarField bubbleRadius = max(bubbleRadiusMin,bubbleRadiusMax*sigma[ip]*sigma[jp]);
            drag[ip] -= 0.5*fvc::interpolate(sigma[jp])*dragCoefficient/fvc::interpolate(bubbleRadius) * mag(Uf[ip]-Uf[jp]) * (Uf[ip]-Uf[jp]);
        }
    }
}
