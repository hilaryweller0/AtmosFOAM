#include "divergentVelocityField.H"

void divergentVelocityField::applyToInternalField(surfaceScalarField& phi) const
{
    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        phi[faceI] = velocityAt(mesh.Cf()[faceI], phi.time()) & mesh.Sf()[faceI];
    }
}

void divergentVelocityField::applyToBoundary(surfaceScalarField& phi, const label patchI) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        bf[faceI] = velocityAt(mesh.Cf()[faceI], phi.time()) & mesh.Sf()[faceI];
    }
}

