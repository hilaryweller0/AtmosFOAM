// Update the velocity in each partition
for (int ucorr = 0; ucorr < nUCorr; ucorr++)
{
    surfaceScalarField ggradh = g*fvc::snGrad(hSum)*mesh.magSf();
    //surfaceScalarField ggradh = lengthScale*g*fvc::snGrad(hSum)/fvc::interpolate(hSum)*mesh.magSf();
    
    //Update prognostic variables.
    for(label ip = 0; ip < nParts; ip++)
    {
        #include "drag.H"
        
        volFlux[ip] = volFlux[ip].oldTime() - dt*
        (
            ((Uf[ip]&fvc::interpolate(fvc::grad(Uf[ip])))&mesh.Sf())
          + ggradh
        );
        
        if (useCoriolis) { volFlux[ip] += -dt*((twoOmegaf^Uf[ip]) & mesh.Sf()); }
        if (useDrag)     { volFlux[ip] +=  dt*(drag[ip] & mesh.Sf()); }
        if (useBuoyancy) { volFlux[ip] +=  dt*buoyancyMagnitude*fvc::interpolate(momentumSource[ip])*(yNorm & mesh.Sf()); }
        if (useGravity)  { volFlux[ip] += -dt*fvc::interpolate(gravity)*(yNorm & mesh.Sf()); }
        
        volScalarField total_source = 0*h[ip];
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                total_source = 0*h[ip];
                if (useMassTransfer) { total_source += dt*sourceMag*sink[jp]*h[jp].oldTime(); }
                if (useDiffusion)    { total_source += max(dt*diffusionCoeff*fvc::laplacian(h[ip].oldTime() - h[jp].oldTime()),0*h[ip]); }
                
                total_source = min(total_source, h[jp].oldTime());
                
                volFlux[ip] += fvc::interpolate(total_source/hSum)*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
            }
        }
        
        /*if (useDiffusion)
        {
            for (label jp = 0; jp < nParts; jp++)
            {
                if (jp != ip)
                {

                    volFlux[ip] += dt*diffusionCoeff*
                               (
                                   ( 

                                       fvc::interpolate(sigma[jp])/fvc::interpolate(sigma[ip])
                                       *fvc::interpolate(fvc::laplacian(sigma[jp]))*(Uf[jp]-Uf[ip])
                             //- fvc::interpolate(fvc::laplacian(sigma[ip].oldTime()))*Uf[ip]
                                   )
                                   & mesh.Sf()
                               );
                }
            }

        }
        
        if (useMassTransfer)
        {
            for (label jp = 0; jp < nParts; jp++)
            {
                if (jp != ip)
                {
                    volFlux[ip] += dt*sourceMag*fvc::interpolate(sink[jp])*
                                      fvc::interpolate(sigma[jp])/fvc::interpolate(sigma[ip])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
                }
            }
        }*/
        
    
        
    
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }
}
