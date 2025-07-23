#include "divergentVelocityField.H"
#include "fvScalarMatrix.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"

namespace Foam
{
divergentVelocityField::divergentVelocityField
(
    const dictionary& dict,
    const bool applyProjection__
)
:
    applyProjection(applyProjection__),
    endTime_(dict.lookupOrDefault<dimensionedScalar>("endTime", scalar(0)))
    {}

void divergentVelocityField::applyToInternalField
(
    surfaceScalarField& phi,
    scalar time
) const
{
    phi = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        phi[faceI] = velocityAt(mesh.Cf()[faceI], time) & mesh.Sf()[faceI];
    }
}

void divergentVelocityField::applyToBoundary
(
    surfaceScalarField& phi,
    const label patchI,
    scalar time
) const
{
    const fvPatch& pat = phi.mesh().boundary()[patchI];
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        bf[faceI] = velocityAt(pat.Cf()[faceI], time) & pat.Sf()[faceI];
    }
}

void divergentVelocityField::project(surfaceScalarField& phi) const
{
    if (applyProjection)
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
