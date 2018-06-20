//Apply mass transfers using operator splitting.
//CURRENTLY ONLY VALID FOR 2 PARTITIONS
if ((useMassTransfer) or (useDiffusion))
{
    for(label ip = 0; ip < nParts; ip++)
    {
        totalSink[ip] *= 0;
        sigmah_temp[ip] = sigmah[ip];
        flux_temp[ip] = flux[ip];
        volFlux_temp[ip] = volFlux[ip];
        hu[ip] = sigmah_temp[ip]*u[ip];
        hu_temp[ip] = sigmah_temp[ip]*u[ip];
        u_temp[ip] = u[ip];
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
                        diffusionCoeff*fvc::laplacian(sigmah_temp[jp] - sigmah_temp[ip]), 
                        minimumMassTransfer
                    ) / (sigmah_temp[ip] + delta); 
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
        if (centreMomentumTransfer)
        {
            sigmaf.updateSum();
            sigma.updateSum();
            sigmah.updateSum();
            volVectorField hu1("hu1", hu.updateSum());
            hu1.write();
            volVectorField u1("u1", u.updateSum());
            u1.write();
            surfaceScalarField volFlux1("volFlux1", volFlux.updateSum());
            volFlux1.write();
            surfaceScalarField flux1("flux1", flux.updateSum());
            flux1.write();

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
                Uf[ip] = fvc::interpolate(u[ip]);
                //volFlux[ip] = Uf[ip] & mesh.Sf();
                //flux[ip] = fvc::interpolate(hu[ip]) & mesh.Sf();
                
                flux[ip] = fvc::interpolate(hu[ip]) & mesh.Sf();
                volFlux[ip] = flux[ip] / (fvc::interpolate(sigmah[ip]) + delta);
                
                //Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                //          *mesh.Sf()/sqr(mesh.magSf());
            }
            volVectorField hu2("hu2", hu.updateSum());
            hu2.write();
            volVectorField u2("u2", u.updateSum());
            u2.write();
            surfaceScalarField volFlux2("volFlux2", volFlux.updateSum());
            volFlux2.write();
            surfaceScalarField flux2("flux2", flux.updateSum());
            flux2.write();
        }
        else if (faceMomentumTransfer)
        {
            for(label ip = 0; ip < nParts; ip++)
            {
                flux[ip] -= dt*fvc::interpolate(totalSink[ip])*flux_temp[ip];
                for (label jp = 0; jp < nParts; jp++)
                {
                    if (jp != ip)
                    {
                        flux[ip] += dt*fvc::interpolate(totalSink[jp])*flux_temp[jp];
                    }
                }
                volFlux[ip] = flux[ip]/fvc::interpolate(sigmah[ip] + delta);
                u[ip] = fvc::reconstruct(volFlux[ip]);
                Uf[ip] = fvc::interpolate(u[ip]);
                Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                          *mesh.Sf()/sqr(mesh.magSf());
            }
        }
    }

    #include "energyTransfers2.H"
}
