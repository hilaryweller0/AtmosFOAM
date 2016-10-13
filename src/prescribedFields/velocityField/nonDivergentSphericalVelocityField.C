#include "nonDivergentSphericalVelocityField.H"
#include "polarPoint.H"

void nonDivergentSphericalVelocityField::applyToInternalField(surfaceScalarField& phi) const
{
    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, mesh, phi.time());
    }
}

void nonDivergentSphericalVelocityField::applyToBoundary(surfaceScalarField& phi, const label patchI) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        bf[faceI] = faceFlux(f, mesh, phi.time());
    }
}

scalar nonDivergentSphericalVelocityField::faceFlux(const face& f, const fvMesh& mesh, const Time& t) const
{
    point p0 = mesh.points()[f.last()];
    point p1 = p0;

    scalar flux = 0;

    scalar minR = -1;
    scalar maxR = -1;
    forAll(f, ip)
    {
        p0 = p1;
        p1 = mesh.points()[f[ip]];
        point pmid = 0.5*(p0 + p1);
        
        if (minR < 0 || mag(p0) < minR) minR = mag(p0);
        if (mag(p0) > maxR) maxR = mag(p0);

        flux += streamfunctionAt(pmid, t) & ((p0 - p1)/mag(p0-p1));
    }

    flux *= (sqr(maxR) - sqr(minR))/2;

    return flux;
}
