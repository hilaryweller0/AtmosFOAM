//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    #include "energyTransfers1.H"
    
    for(label ip = 0; ip < nParts; ip++)
    {
        sigmah_temp[ip] = sigmah[ip];
    }
    
    for(label ip = 0; ip < nParts; ip++)
    {
        if (useMassTransfer) 
        { 
            totalSink[ip] = min( sourceMag*sink[ip], 1./dt );
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
        sigmah[ip] -= dt*totalSink[ip]*sigmah_temp[ip];
        
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                sigmah[ip] += dt*totalSink[jp]*sigmah_temp[jp];
            }
        }
    }
    //Update total height.
    //sigmah.updateSum();
    h = sigmah.sum();

    //Add momentum transfer from mass transfer terms.
    for(label ip = 0; ip < nParts; ip++)
    {
        flux_temp[ip] = flux[ip];
    }
    
    //volFlux_temp = volFlux;
    for(label ip = 0; ip < nParts; ip++)
    {
        flux[ip] -= dt*fvc::interpolate(totalSink[ip])*flux_temp[ip];
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                //volFlux[ip] += fvc::interpolate(totalSink[jp]/sigmah_temp[jp])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
                flux[ip] += dt*fvc::interpolate(totalSink[jp])*flux_temp[jp];
            }
        }
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] = flux[ip]/fvc::interpolate(sigmah[ip]);
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        //Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
        //        *mesh.Sf()/sqr(mesh.magSf());
    }

    #include "energyTransfers2.H"
}
