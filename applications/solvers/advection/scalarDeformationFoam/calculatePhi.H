// Update phi (and U and Uf) based on runTime.time()
// phi assumes the domain in a cylinder rotated about the z axis

{
    const scalar t = runTime.time().value();
    phi.internalField() = 0;

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
            const scalar r = 0.5*
            (
                Foam::sqrt(sqr(p0.x()) + sqr(p0.y()))
              + Foam::sqrt(sqr(p1.x()) + sqr(p1.y()))
            );
            point pmid = 0.5*(p0 + p1);
            scalar z = pmid.z();
            vector normal = pmid - z*vector(0,0,1);
            normal /= mag(normal);
            pmid /= mag(pmid);
            pmid = z*vector(0,0,1) + r*normal;
            scalar x = r*Foam::atan2(pmid.y(), pmid.x());
            vector streamFunction = normal*
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
