// Update phi (and U and Uf) based on runTime.time() on a plane
// for deformational flow 

{
    const scalar t = runTime.time().value();
    phi.ref() = 0;

    // Set all internal values of phi
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        point p0 = mesh.points()[f.last()];
        point p1 = p0;
        
        // Loop around all edges of the face to set the flux from the 
        // streamfunction using Stokes circulation theorem
        forAll(f, ip)
        {
            p0 = p1;
            p1 = mesh.points()[f[ip]];
            point pmid = 0.5*(p0 + p1);
            scalar z = pmid.z();
            scalar x = pmid.x();
            vector streamFunction = vector(0,-1,0)*
            (
                psiScale*sqr(0.5*L.x()/M_PI)/endTime*sqr
                (
                    Foam::sin(2*M_PI*(x/L.x() - t/endTime))
                   *Foam::cos(M_PI*z/L.z())
                )
                *Foam::cos(M_PI*t/endTime)
              - L.x()*z/endTime
            );
            phi[faceI] += streamFunction & (p0 - p1);
        }
    }
    
    // Update U and Uf for post processing
    U = fvc::reconstruct(phi);
    Uf = linearInterpolate(U);
    Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
}
