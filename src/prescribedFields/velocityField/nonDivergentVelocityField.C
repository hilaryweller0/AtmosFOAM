#include "nonDivergentVelocityField.H"

namespace Foam
{

nonDivergentVelocityField::nonDivergentVelocityField(const dictionary& dict)
:
    endTime_(dict.lookupOrDefault<dimensionedScalar>("endTime", scalar(0)))
{}


void nonDivergentVelocityField::applyToInternalField
(
    surfaceScalarField& phi,
    scalar time
) const
{
    phi = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, mesh, time);
    }
}

void nonDivergentVelocityField::applyToBoundary
(
    surfaceScalarField& phi,
    const label patchI,
    scalar time
) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        bf[faceI] = faceFlux(f, mesh, time);
    }
}

scalar nonDivergentVelocityField::faceFlux
(
    const face& f,
    const fvMesh& mesh,
    scalar time
) const
{
    point p0 = mesh.points()[f.last()];
    point p1;
    vector s0 = streamfunctionAt(p0, time);
    vector s1;

    scalar flux = 0;

    forAll(f, ip)
    {
        p1 = mesh.points()[f[ip]];
        s1 = streamfunctionAt(p1, time);
        //point pmid = 0.5*(p0 + p1);
        //flux += streamfunctionAt(pmid, t) & (p0 - p1);
        flux += (s0 + s1) & (p0 - p1);
        p0 = p1;
        s0 = s1;
    }

    return 0.5*flux;
}

void nonDivergentVelocityField::project(surfaceScalarField& phi) const
{}

}
