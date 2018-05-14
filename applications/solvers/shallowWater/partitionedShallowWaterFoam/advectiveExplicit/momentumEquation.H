// Update the velocity in each partition
for (int ucorr = 0; ucorr < nUCorr; ucorr++)
{
    surfaceScalarField ggradh = g*fvc::snGrad(h)*mesh.magSf();
    
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
          + (1-offCentre)*g*fvc::snGrad(h.oldTime())*mesh.magSf()
        );
        
        if (useCoriolis) { volFlux[ip] += -dt*((twoOmegaf^Uf[ip]) & mesh.Sf()); }
        if (useDrag)     { volFlux[ip] +=  dt*(drag[ip] & mesh.Sf()); }
        if (useBuoyancy) 
        { 
            volFlux[ip] += dt*
            (
                //Foam::sin(M_PI*runTime.time().value()/100.)*Foam::sin(M_PI*runTime.time().value()/100.)*
                buoyancyMagnitude*fvc::interpolate(momentumSource[ip])*(yNorm & mesh.Sf())
            ); 
        }
        if (useGravity)  { volFlux[ip] += -dt*fvc::interpolate(gravity)*(yNorm & mesh.Sf()); }
        
        flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
        
        /*u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());*/
    }
    flux.updateSum();
    //Helmholtz equation and final update of height gradient term.
    if (semiImplicit)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(h)
          + fvc::div(flux.sum())
          - fvm::laplacian(offCentre*dt*g*h.oldTime(), h)
        );
        hEqn.solve();
    }
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] += -offCentre*dt*g*fvc::snGrad(h)*mesh.magSf();
    }
    
    for(label ip = 0; ip < nParts; ip++)
    {
        flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                *mesh.Sf()/sqr(mesh.magSf());
    }
}
