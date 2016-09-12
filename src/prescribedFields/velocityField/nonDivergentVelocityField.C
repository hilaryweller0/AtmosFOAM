#include "nonDivergentVelocityField.H"

void nonDivergentVelocityField::applyToInternalField(surfaceScalarField& phi) const
{
    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, mesh, phi.time());
    }
}

void nonDivergentVelocityField::applyToBoundary(surfaceScalarField& phi, const label patchI) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        bf[faceI] = faceFlux(f, mesh, phi.time());
    }
}

scalar nonDivergentVelocityField::faceFlux(const face& f, const fvMesh& mesh, const Time& t) const
{
    point p0 = mesh.points()[f.last()];
    point p1 = p0;

    scalar flux = 0;

    forAll(f, ip)
    {
        p0 = p1;
        p1 = mesh.points()[f[ip]];
        point pmid = 0.5*(p0 + p1);
        flux += streamfunctionAt(pmid, t) & (p0 - p1);
    }

    return flux;
}
