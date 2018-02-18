// Update the velocity in each partition
for (int ucorr = 0; ucorr < nUCorr; ucorr++)
{
    surfaceScalarField ggradh = g*fvc::snGrad(hSum)*mesh.magSf();
    //surfaceScalarField ggradh = lengthScale*g*fvc::snGrad(hSum)/fvc::interpolate(hSum)*mesh.magSf();
    
    fluxSum *= 0;
    //Update prognostic variables.
    for(label ip = 0; ip < nParts; ip++)
    {
        #include "drag.H"
        
        volFlux[ip] = volFlux[ip].oldTime() - dt*
        (
            (1-offCentre) *
            (
                ( Uf[ip].oldTime() & fvc::interpolate( fvc::grad(Uf[ip].oldTime()) ) ) & mesh.Sf()
            )
          + offCentre *
            (
                ( Uf[ip] & fvc::interpolate( fvc::grad(Uf[ip]) ) ) & mesh.Sf()
            )
          + (1-offCentre)*g*fvc::snGrad(hSum.oldTime())*mesh.magSf()
        );
        
        if (useCoriolis) { volFlux[ip] += -dt*((twoOmegaf^Uf[ip]) & mesh.Sf()); }
        if (useDrag)     { volFlux[ip] +=  dt*(drag[ip] & mesh.Sf()); }
        if (useBuoyancy) { volFlux[ip] +=  dt*buoyancyMagnitude*fvc::interpolate(momentumSource[ip])*(yNorm & mesh.Sf()); }
        if (useGravity)  { volFlux[ip] += -dt*fvc::interpolate(gravity)*(yNorm & mesh.Sf()); }
        
        fluxSum += fvc::interpolate(h[ip])*volFlux[ip];
        
        /*u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());*/
    }
    
    //Helmholtz equation and final update of height gradient term.
    if (semiImplicit)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(hSum)
          + fvc::div(fluxSum)
          - fvm::laplacian(offCentre*dt*g*hSum.oldTime(), hSum)
        );
        hEqn.solve();
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] += -offCentre*dt*g*fvc::snGrad(hSum)*mesh.magSf();
    }
    
    //Add momentum transfer from mass transfer terms.
    for(label ip = 0; ip < nParts; ip++)
    {
        for (label jp = 0; jp < nParts; jp++)
        {
            if (jp != ip)
            {
                volFlux[ip] += fvc::interpolate(total_sink[jp]/h_temp[jp])*( (Uf[jp]-Uf[ip]) & mesh.Sf() );
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
