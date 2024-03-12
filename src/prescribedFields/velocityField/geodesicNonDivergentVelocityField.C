#include "geodesicNonDivergentVelocityField.H"

namespace Foam
{
geodesicNonDivergentVelocityField::geodesicNonDivergentVelocityField
(
    const dictionary& dict
)
:
    earthRadius_(readScalar(dict.lookup("earthRadius")))
{}

void geodesicNonDivergentVelocityField::applyToInternalField
(
    surfaceScalarField& phi
) const
{
    const fvMesh& mesh = phi.mesh();

    // Get reference to spherical geometry
    const sphericalMeshData& spherical = sphericalMeshData::New
    (
        mesh, earthRadius_
    );

    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, spherical, phi.time());
    }
}

void geodesicNonDivergentVelocityField::applyToBoundary
(
    surfaceScalarField& phi, const label patchI
) const
{
    const fvMesh& mesh = phi.mesh();
    // Get reference to spherical geometry
    const sphericalMeshData& spherical = sphericalMeshData::New
    (
        mesh, earthRadius_
    );

    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        bf[faceI] = faceFlux(f, spherical, phi.time());
    }
}

scalar geodesicNonDivergentVelocityField::faceFlux
(
    const face& f,
    const sphericalMeshData& spherical,
    const Time& t
) const
{
    point p0 = spherical.mesh().points()[f.last()];
    point p1;
    vector s0 = streamfunctionAt(f.last(), spherical, t);
    vector s1;
    scalar mags0 = mag(s0);
    scalar mags1;

    scalar flux = 0;

    forAll(f, ip)
    {
        p1 = spherical.mesh().points()[f[ip]];
        s1 = streamfunctionAt(f[ip], spherical, t);
        mags1 = mag(s1);
        
        // Average streamfunction is the average direction and the average
        // magnitude
        vector s = s0/mags0 + s1/mags1;
        s *= 0.5*(mags0 + mags1)/mag(s);
        
        flux += s & (p0 - p1);
        p0 = p1;
        s0 = s1;
        mags0 = mags1;
    }

/*    scalar r0 = mag(spherical.mesh().points()[f[0]]);
    scalar r1 = mag(spherical.mesh().points()[f[1]]);
    scalar r2 = mag(spherical.mesh().points()[f[2]]);
    scalar r3 = mag(spherical.mesh().points()[f[3]]);
    if (r0 < 1 && r1 < 1 && r2 < 1 && r3 < 1 && mag(flux) > 1e-6)
    {
        Info << "Flux = " << flux << endl;
        forAll(f, ip)
        {
            p1 = spherical.mesh().points()[f[ip]];
            s1 = streamfunctionAt(f[ip], spherical, t);
            mags1 = mag(s1);
            
            // Average streamfunction is the average direction and the average
            // magnitude
            vector s = s0/mags0 + s1/mags1;
            s *= 0.5*(mags0 + mags1)/mag(s);
            
            vector pmidHat = 0.5*(p0 + p1);
            pmidHat /= mag(pmidHat);
            
            Info << "pmidHat= " << pmidHat
                 << " sHat= " << s/mag(s) 
                 << " s0 = " << s0
                 << " s1 = " << s1 << nl
                 << "p0hat = " << p0/mag(p0) << " p1hat = " << p1/mag(p1)
                 << " s&(p0-p1)= " << (s & (p0 - p1)) << endl;
                
            p0 = p1;
            s0 = s1;
            mags0 = mags1;
        }
    }
*/
    return flux;
}

void geodesicNonDivergentVelocityField::project(surfaceScalarField& phi) const
{}
}
