// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    h[ip] = h[ip].oldTime() - dt*
    (
        fvc::div(volFlux[ip], h[ip])
    );
    
    //Calculate total height.
    if (ip == 0) hSum  = h[ip];
    else         hSum += h[ip];
}

//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    for(label ip = 0; ip < nParts; ip++)
    {
        total_sink[ip]   = 0*h[ip];
        //volScalarField total_source = 0*h[ip];
        if (useMassTransfer) { total_sink[ip] += dt*sourceMag*sink[ip]*h[ip]; }
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                if (useDiffusion)    { total_sink[ip] += max(dt*diffusionCoeff*fvc::laplacian(h[jp] - h[ip]),0*h[ip]); }
                
                //total_source = 0*h[ip];
                //if (useMassTransfer) { total_source += dt*sourceMag*sink[jp]*h[jp].oldTime(); }
                //if (useDiffusion)    { total_source += max(dt*diffusionCoeff*fvc::laplacian(h[ip].oldTime() - h[jp].oldTime()),0*h[ip]); }
                
                //total_source = min(total_source, h[jp].oldTime());
                
                //h[ip] += total_source;
            }
        }
        //h[ip] += -total_sink;
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        total_sink[ip] = min(total_sink[ip], h[ip]);
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        h[ip] -= total_sink[ip];
        
        for (label jp = 0; jp < nParts; jp++)
        {
            h[ip] += total_sink[jp];
        }
        
        //Calculate total height.
        if (ip == 0) hSum  = h[ip];
        else         hSum += h[ip];
    }
}
