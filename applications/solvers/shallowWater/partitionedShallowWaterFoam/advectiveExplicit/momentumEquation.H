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
        
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }
    
    for(label ip = 0; ip < nParts; ip++)
    {
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                volFlux[ip] += fvc::interpolate(total_sink[ip]/h[jp])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
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
}
