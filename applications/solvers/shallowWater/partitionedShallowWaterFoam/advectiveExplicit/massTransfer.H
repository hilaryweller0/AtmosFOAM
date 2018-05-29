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
            totalSink[ip] += sourceMag*sink[ip];
            //totalSink[ip] = min( sourceMag*sink[ip], 1./dt );
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
    h = sigmah.sum();

    //Add momentum transfer from mass transfer terms.
    for(label ip = 0; ip < nParts; ip++)
    {
        flux_temp[ip] = flux[ip];
        volFlux_temp[ip] = volFlux[ip];
    }
    
    if (implicitMomentumTransfer)
    {
        flux[0] = flux_temp[0] * (1 + dt*fvc::interpolate(totalSink[1])) / 
                    (1 + dt*(fvc::interpolate(totalSink[0]) + fvc::interpolate(totalSink[1])))
                  + flux_temp[1] * dt*fvc::interpolate(totalSink[1]) / 
                    (1 + dt*(fvc::interpolate(totalSink[0]) + fvc::interpolate(totalSink[1])));
        flux[1] = flux_temp[1] * (1 + dt*fvc::interpolate(totalSink[0])) / 
                    (1 + dt*(fvc::interpolate(totalSink[0]) + fvc::interpolate(totalSink[1])))
                  + flux_temp[0] * dt*fvc::interpolate(totalSink[0]) / 
                    (1 + dt*(fvc::interpolate(totalSink[0]) + fvc::interpolate(totalSink[1])));
    }
    else
    {
        surfaceScalarField denominator("denominator",1 + dt*
            (
                fvc::interpolate(totalSink[0])*fvc::interpolate(sigmah_temp[0])/fvc::interpolate(sigmah_temp[1] + delta)
              + fvc::interpolate(totalSink[1])*fvc::interpolate(sigmah_temp[1])/fvc::interpolate(sigmah_temp[0] + delta)
            ));
        volFlux[0] = volFlux_temp[0] * (1 + dt*fvc::interpolate(totalSink[0])*fvc::interpolate(sigmah_temp[0])/fvc::interpolate(sigmah_temp[1] + delta)) / denominator
        + volFlux_temp[1] * dt*fvc::interpolate(totalSink[1])*fvc::interpolate(sigmah_temp[1])/fvc::interpolate(sigmah_temp[0] + delta) / denominator;
        volFlux[1] = volFlux_temp[1] * (1 + dt*fvc::interpolate(totalSink[1])*fvc::interpolate(sigmah_temp[1])/fvc::interpolate(sigmah_temp[0] + delta)) / denominator
        + volFlux_temp[0] * fvc::interpolate(totalSink[0])*fvc::interpolate(sigmah_temp[0])/fvc::interpolate(sigmah_temp[1] + delta) / denominator;
    }   
    /*else
    {
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
    }*/

    for(label ip = 0; ip < nParts; ip++)
    {
        //volFlux[ip] = flux[ip]/fvc::interpolate(sigmah[ip]);
        flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        //Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
        //        *mesh.Sf()/sqr(mesh.magSf());
    }

    #include "energyTransfers2.H"
}
