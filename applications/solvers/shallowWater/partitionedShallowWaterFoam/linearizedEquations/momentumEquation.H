// Update the velocity in each partition
for (int ucorr = 0; ucorr < nUCorr; ucorr++)
{
    //Update prognostic variables.
    for(label ip = 0; ip < nParts; ip++)
    {
        volFlux[ip] = volFlux[ip].oldTime() - dt*
        (
            ( ( UF[ip] & fvc::interpolate( fvc::grad(Uf[ip]) ) ) & mesh.Sf() )
          + gravity*fvc::snGrad(h)*mesh.magSf()
        );
        
        /*Uf[ip] = Uf[ip].oldTime() - dt*
        (
            ( UF[ip] & fvc::interpolate( fvc::grad(Uf[ip]) ) )
          + gravity*fvc::snGrad(h) * mesh.Sf() / mesh.magSf()
        );*/
    }
    
    for(label ip = 0; ip < nParts; ip++)
    {
        //volFlux[ip] = Uf[ip] & mesh.Sf();
        u[ip] = fvc::reconstruct(volFlux[ip]);
        Uf[ip] = fvc::interpolate(u[ip]);
        
        //Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
        //        *mesh.Sf()/sqr(mesh.magSf());
    }
}
