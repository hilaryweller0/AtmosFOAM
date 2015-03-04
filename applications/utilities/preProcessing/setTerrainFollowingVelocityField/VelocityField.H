#ifndef VELOCITY_FIELD_H
#define VELOCITY_FIELD_H

#include "VelocityProfile.H"
#include "fvCFD.H"

class VelocityField
{
    public:
    VelocityField(const VelocityProfile& profile);

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyTo(GeometricField<Type, PatchField, GeoMesh>& field);

    private:
    const VelocityProfile& profile;

    template<class Type>
    void applyToInternalField(GeometricField<Type, fvsPatchField, surfaceMesh>& field);

    template<class Type>
    void applyToInternalField(GeometricField<Type, fvPatchField, volMesh>& field);

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyToBoundary
    (
        const word name, GeometricField<Type, PatchField, GeoMesh>& field
    );

    label findBoundaryPatchIndex(const fvMesh& mesh, const word& name);
};

#endif
