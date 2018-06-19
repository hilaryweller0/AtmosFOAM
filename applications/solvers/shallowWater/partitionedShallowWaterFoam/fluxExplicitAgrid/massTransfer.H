//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    if (runTime.time().value() > 100)
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            sink[0] *= 0;
            momentumSource[ip] *= 0;
        }
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        totalSink[ip] *= 0;
        sigmah_temp[ip] = sigmah[ip];
        hu[ip] = sigmah_temp[ip]*u[ip];
        hu_temp[ip] = sigmah_temp[ip]*u[ip];
    }
    
    #include "energyTransfers1.H"
    
    for(label ip = 0; ip < nParts; ip++)
    {
        if (useMassTransfer) 
        { 
            totalSink[ip] += sourceMag*sink[ip];
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
                    ) / (sigmah_temp[ip] + scalar(1e-16)); 
                }
            }
        }
        
        if ( !(implicitMassTransfer) or !(implicitMomentumTransfer) )
        {
            totalSink[ip] = min( totalSink[ip], 1./dt );
        }
    }
    
    
    if (implicitMassTransfer)
    {
        sigmah[0] = sigmah_temp[0] * (1 + dt*totalSink[1]) / (1 + dt*(totalSink[0] + totalSink[1]))
                  + sigmah_temp[1] * dt*totalSink[1] / (1 + dt*(totalSink[0] + totalSink[1]));
        sigmah[1] = sigmah_temp[1] * (1 + dt*totalSink[0]) / (1 + dt*(totalSink[0] + totalSink[1]))
                  + sigmah_temp[0] * dt*totalSink[0] / (1 + dt*(totalSink[0] + totalSink[1]));
    }
    else
    {
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
    }
    //Update total height.
    //sigmah.updateSum();
    h = sigmah[0] + sigmah[1];

    //Add momentum transfer from mass transfer terms.
    if (implicitMomentumTransfer)
    {
        hu[0] = hu_temp[0] * (1 + dt*totalSink[1]) / 
                    (1 + dt*(totalSink[0] + totalSink[1]))
                  + hu_temp[1] * dt*totalSink[1] / 
                    (1 + dt*(totalSink[0] + totalSink[1]));
        hu[1] = hu_temp[1] * (1 + dt*totalSink[0]) / 
                    (1 + dt*(totalSink[0] + totalSink[1]))
                  + hu_temp[0] * dt*totalSink[0] / 
                    (1 + dt*(totalSink[0] + totalSink[1]));
        u[0] = hu[0] / (sigmah[0] + delta);
        u[1] = hu[1] / (sigmah[1] + delta);
    }
    else
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            hu[ip] -= dt*totalSink[ip]*hu_temp[ip];
            for (label jp = 0; jp < nParts; jp++)
            {
                if (jp != ip)
                {
                    hu[ip] += dt*totalSink[jp]*hu_temp[jp];
                }
            }
            u[ip] = hu[ip] / (sigmah[ip] + delta);
        }
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        Uf[ip] = fvc::interpolate(u[ip]);
        
        volFlux[ip] = Uf[ip] & mesh.Sf();
        flux[ip] = fvc::interpolate(hu[ip]) & mesh.Sf();
        
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }

    #include "energyTransfers2.H"
}
