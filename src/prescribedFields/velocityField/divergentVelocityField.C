#include "divergentVelocityField.H"

divergentVelocityField::divergentVelocityField(const bool divFree__)
:
    divFree(divFree__)
{}

void divergentVelocityField::applyTo(surfaceScalarField& phi) const
{
    velocityField::applyTo(phi);
    if (divFree)
    {
        project(phi);
    }
    else Info << "Not projecting" << endl;
}

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
    const fvPatch& pat = phi.mesh().boundary()[patchI];
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        bf[faceI] = velocityAt(pat.Cf()[faceI], phi.time()) & pat.Sf()[faceI];
    }
}

void divergentVelocityField::project(surfaceScalarField& phi) const
{
    Info << "Projecting" << endl;

}
