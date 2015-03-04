#include "VelocityField.H"

VelocityField::VelocityField(const VelocityProfile& profile) : profile(profile) {};

template<class Type, template<class> class PatchField, class GeoMesh>
void VelocityField::applyTo(GeometricField<Type, PatchField, GeoMesh>& field)
{
    applyToInternalField(field);
    applyToBoundary("inlet", field);
    applyToBoundary("outlet", field);
}

template<class Type>
void VelocityField::applyToInternalField(GeometricField<Type, fvsPatchField, surfaceMesh>& field)
{
    forAll(field, faceI)
    {
        const point& Cf = field.mesh().Cf()[faceI];
        field[faceI] = profile.velocityAt(Cf);
    }
}

template<class Type>
void VelocityField::applyToInternalField(GeometricField<Type, fvPatchField, volMesh>& field)
{
    forAll(field, cellI)
    {
        const point& C = field.mesh().C()[cellI];
        field[cellI] = profile.velocityAt(C);
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void VelocityField::applyToBoundary
(
    const word name, GeometricField<Type, PatchField, GeoMesh>& field
)
{
    label boundaryI = findBoundaryPatchIndex(field.mesh(), name);
    forAll(field.boundaryField()[boundaryI], cellI)
    {
        const point& face = field.mesh().Cf().boundaryField()[boundaryI][cellI];
        field.boundaryField()[boundaryI][cellI] = profile.velocityAt(face);
    }
}

label VelocityField::findBoundaryPatchIndex(const fvMesh& mesh, const word& name)
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].name() == name)
        {
            return patchI;
        }
    }
    FatalErrorIn("VelocityField")
        << " no boundary called " << name << ". The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);

    return -1;
}
