if (runTime.time().value() > 100)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        sink[ip] *= 0;
        momentumSink[ip] *= 0;
        momentumSource[ip] *= 0;
    }
}

// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    sigmah_temp[ip] = sigmah[ip].oldTime();
    sigmah[ip] = sigmah_temp[ip] - dt*
    (
        (1-offCentre)*fvc::div(volFlux[ip].oldTime(), sigmah[ip].oldTime())
      + offCentre*fvc::div(volFlux[ip], sigmah[ip])
    );
    
    if (useOperatorSplitting)
    {
        sigmah_temp[ip] = sigmah[ip];
    }
}
//Update total height.
sigmah.updateSum();
h = sigmah.sum();

//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    for(label ip = 0; ip < nParts; ip++)
    {
        totalSink[ip]   = 0*sigmah[ip];
        if (useMassTransfer) 
        { 
            totalSink[ip] += dt*sourceMag*sink[ip]*sigmah_temp[ip]; 
        }
        
        if (useDiffusion)
        {
            for (label jp = 0; jp < nParts; jp++)
            {
                if (jp != ip)
                {
                    totalSink[ip] += max
                    (
                        dt*diffusionCoeff*fvc::laplacian(sigmah_temp[jp] - sigmah_temp[ip]), 
                        minimumMassTransfer
                    ); 
                }
            }
        }
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        totalSink[ip] = min(totalSink[ip], sigmah_temp[ip]);
    }
    
    for (label ip = 0; ip < nParts; ip++)
    {
        sigmah[ip] -= totalSink[ip];
        
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                sigmah[ip] += totalSink[jp];
            }
        }
    }
    //Update total height.
    //sigmah.updateSum();
    h = sigmah.sum();
}
