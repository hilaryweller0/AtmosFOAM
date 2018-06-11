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
    /*else
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
    }*/
    else
    {
        for(label ip = 0; ip < nParts; ip++)
        {
            //flux[ip] -= dt*fvc::interpolate(totalSink[ip])*flux_temp[ip];
            hu[ip] -= dt*totalSink[ip]*hu_temp[ip];
            for (label jp = 0; jp < nParts; jp++)
            {
                if (jp != ip)
                {
                    //volFlux[ip] += dt*fvc::interpolate(totalSink[jp])*fvc::interpolate(sigmah_temp[jp])/fvc::interpolate(sigmah_temp[ip])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
                    //flux[ip] += dt*fvc::interpolate(totalSink[jp])*flux_temp[jp];
                    hu[ip] += dt*totalSink[jp]*hu_temp[jp];
                }
            }
            //u[ip] = hu[ip] / (sigmah[ip] + delta);
        }
        
        //Correction to keep u bounded for sharp gradients of the transfer term.
        sigmah_temp[0] = sigmah[0];
        /*forAll(u[0],celli)
        {
            uminx = min(u_temp[0][celli].component(vector::X),u_temp[1][celli].component(vector::X));
            umaxx = max(u_temp[0][celli].component(vector::X),u_temp[1][celli].component(vector::X));
            uminy = min(u_temp[0][celli].component(vector::Y),u_temp[1][celli].component(vector::Y));
            umaxy = max(u_temp[0][celli].component(vector::Y),u_temp[1][celli].component(vector::Y));
            uminz = min(u_temp[0][celli].component(vector::Z),u_temp[1][celli].component(vector::Z));
            umaxz = max(u_temp[0][celli].component(vector::Z),u_temp[1][celli].component(vector::Z));
            
            u0x = max( min(u[0][celli].component(vector::X), umaxx), uminx );
            u1x = max( min(u[1][celli].component(vector::X), umaxx), uminx );
            u0y = max( min(u[0][celli].component(vector::Y), umaxy), uminy );
            u1y = max( min(u[1][celli].component(vector::Y), umaxy), uminy );
            u0z = max( min(u[0][celli].component(vector::Z), umaxz), uminz );
            u1z = max( min(u[1][celli].component(vector::Z), umaxz), uminz );
            
            //umin = vector(uminx,uminy,uminz);
            //umax = vector(umaxx,umaxy,umaxz);
            
            //u[0][celli] = max( min(u[0][celli], umax), umin );
            //u[1][celli] = max( min(u[1][celli], umax), umin );
            
            u[0][celli] = vector(u0x,u0y,u0z);
            u[1][celli] = vector(u1x,u1y,u1z);
        
            /*if (mag(u[0][celli] - u[1][celli]) != 0.)
            //{
                sigmah_temp[0][celli] = mag(hu[0][celli] + hu[1][celli] - (sigmah[0][celli] + sigmah[1][celli])*u[1][celli])/mag(u[0][celli] - u[1][celli]);
                sigmah_temp[0][celli] = max( min(sigmah_temp[0][celli], sigmah[0][celli] + sigmah[1][celli] - delta.value()), delta.value());
                if ((sigmah[0][celli] + sigmah[1][celli] - sigmah_temp[0][celli] < 0.) or (sigmah_temp[0][celli] < 0.))
                {
                    Info << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
                    Info << "h0 old/new: " << sigmah[0][celli] << " " << sigmah_temp[0][celli] << endl;
                    Info << "h1 old/new: " << sigmah[1][celli] << " " << sigmah[0][celli] + sigmah[1][celli] - sigmah_temp[0][celli] << endl;
                } 
                sigmah[1][celli] = sigmah[0][celli] + sigmah[1][celli] - sigmah_temp[0][celli];
                sigmah[0][celli] = sigmah_temp[0][celli];
            }
        }*/
        hu[0] = sigmah[0]*u[0];
        hu[1] = sigmah[1]*u[1];
    }

    for(label ip = 0; ip < nParts; ip++)
    {
        //volFlux[ip] = flux[ip]/fvc::interpolate(sigmah[ip] + delta);
        //flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
        //u[ip] = fvc::reconstruct(volFlux[ip]);
        
        
        Uf[ip] = fvc::interpolate(u[ip]);
        
        volFlux[ip] = Uf[ip] & mesh.Sf();
        flux[ip] = fvc::interpolate(hu[ip]) & mesh.Sf();
        
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }

    #include "energyTransfers2.H"
}
