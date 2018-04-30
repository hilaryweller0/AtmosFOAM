//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    #include "energyTransfers1.H"

    for(label ip = 0; ip < nParts; ip++)
    {
        sigmah_temp[ip] = sigmah[ip];
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

    //Add momentum transfer from mass transfer terms.
    for(label ip = 0; ip < nParts; ip++)
    {
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                volFlux[ip] += fvc::interpolate(totalSink[jp]/sigmah_temp[jp])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
            }
        }
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }

    #include "energyTransfers2.H"
}
