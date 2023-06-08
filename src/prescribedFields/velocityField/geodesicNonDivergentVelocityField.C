#include "geodesicNonDivergentVelocityField.H"
#include "fvCFD.H"

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
    point p1 = p0;
    vector s0 = streamfunctionAt(f.last(), spherical, t);
    vector s1 = s0;

    scalar flux = 0;

    forAll(f, ip)
    {
        p0 = p1;
        s0 = s1;
        p1 = spherical.mesh().points()[f[ip]];
        s1 = streamfunctionAt(f[ip], spherical, t);
        flux += (s0 + s1) & (p0 - p1);
    }

    return 0.5*flux;
}

void geodesicNonDivergentVelocityField::project(surfaceScalarField& phi) const
{}
