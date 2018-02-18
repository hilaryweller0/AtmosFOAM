// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    h_temp[ip] = h[ip].oldTime();
    h[ip] = h_temp[ip] - dt*
    (
        (1-offCentre)*fvc::div(volFlux[ip].oldTime(), h[ip].oldTime())
      + offCentre*
        fvc::div(volFlux[ip], h[ip])
    );
    
    if (useOperatorSplitting)
    {
        h_temp[ip] = h[ip];
    }
    
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
        if (useMassTransfer) 
        { 
            total_sink[ip] += dt*sourceMag*sink[ip]*h_temp[ip]; 
        }
        
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                if (useDiffusion)    
                {
                    total_sink[ip] += max
                    (
                        dt*diffusionCoeff*fvc::laplacian(h_temp[jp] - h_temp[ip]), 
                        minimumMassTransfer
                    ); 
                }
            }
        }
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        total_sink[ip] = min(total_sink[ip], h_temp[ip]);
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        h[ip] -= total_sink[ip];
        
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                h[ip] += total_sink[jp];
            }
        }
        
        //Calculate total height.
        if (ip == 0) hSum  = h[ip];
        else         hSum += h[ip];
    }
}
