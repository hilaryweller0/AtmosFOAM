#ifndef VELOCITY_FIELD_H
#define VELOCITY_FIELD_H

#include "VelocityProfile.H"
#include "fvCFD.H"

class VelocityField
{
    public:
    VelocityField(const VelocityProfile& profile) : profile(profile) {};

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyTo(GeometricField<Type, PatchField, GeoMesh>& field) const
    {
        applyToInternalField(field);
        applyToBoundary("inlet", field);
        applyToBoundary("outlet", field);
    }

    private:
    const VelocityProfile& profile;

    template<class Type>
    void applyToInternalField(GeometricField<Type, fvsPatchField, surfaceMesh>& field) const
    {
        forAll(field, faceI)
        {
            const point& Cf = field.mesh().Cf()[faceI];
            field[faceI] = profile.velocityAt(Cf);
        }
    }

    template<class Type>
    void applyToInternalField(GeometricField<Type, fvPatchField, volMesh>& field) const
    {
        forAll(field, cellI)
        {
            const point& C = field.mesh().C()[cellI];
            field[cellI] = profile.velocityAt(C);
        }
    }

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyToBoundary
    (
        const word name, GeometricField<Type, PatchField, GeoMesh>& field
    ) const
    {
        label boundaryI = findBoundaryPatchIndex(field.mesh(), name);
        forAll(field.boundaryField()[boundaryI], cellI)
        {
            const point& face = field.mesh().Cf().boundaryField()[boundaryI][cellI];
            field.boundaryField()[boundaryI][cellI] = profile.velocityAt(face);
        }
    }

    label findBoundaryPatchIndex(const fvMesh& mesh, const word& name) const
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
};

#endif
