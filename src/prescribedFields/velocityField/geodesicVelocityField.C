#include "geodesicVelocityField.H"
#include "sphericalMeshData.H"
#include "fvScalarMatrix.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"

namespace Foam
{
geodesicVelocityField::geodesicVelocityField(const dictionary& dict)
:
    earthRadius_(dict.lookup("earthRadius")),
    endTime_(dict.lookup("endTime_")),
    applyProjection_(readBool(dict.lookup("applyProjection")))
{}

void geodesicVelocityField::applyToInternalField
(
    surfaceScalarField& phi,
    scalar time
) const
{
    // Get reference to spherical geometry
    const sphericalMeshData& spherical = sphericalMeshData::New
    (
        phi.mesh(), earthRadius_.value()
    );

    phi = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        phi[faceI] = velocityAt
        (
            spherical.faceCentres()[faceI],
            time
        ) & mesh.Sf()[faceI];
    }
}

void geodesicVelocityField::applyToBoundary
(
    surfaceScalarField& phi,
    const label patchI,
    scalar time
) const
{
    // Get reference to spherical geometry
    const sphericalMeshData& spherical = sphericalMeshData::New
    (
        phi.mesh(), earthRadius_.value()
    );

    const fvPatch& pat = phi.mesh().boundary()[patchI];
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        label i = pat.start() + faceI;
        bf[faceI] = velocityAt
        (
            spherical.faceCentres()[i], time
        ) & spherical.faceAreas()[i];
    }
}

void geodesicVelocityField::project(surfaceScalarField& phi) const
{
    if (applyProjection_)
    {
        const fvMesh& mesh = phi.mesh();
        volScalarField P
        (
            IOobject("ProjectionP", phi.time().name(), mesh),
            mesh,
            dimensionedScalar("P", phi.dimensions()/dimLength, scalar(0)),
            "zeroGradient"
        );
        fvScalarMatrix PEqn(fvm::laplacian(P) + fvc::div(phi));
        PEqn.setReference(0,0);
        PEqn.solve();
        phi += PEqn.flux();
    }
}
}
