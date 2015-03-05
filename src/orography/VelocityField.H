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
        applyToBoundary("top", field);
        applyToBoundary("ground", field);
    }

    private:
    const VelocityProfile& profile;

    void applyToInternalField(surfaceVectorField& field) const
    {
        forAll(field, faceI)
        {
            const point& Cf = field.mesh().Cf()[faceI];
            field[faceI] = profile.velocityAt(Cf);
        }
    }

    void applyToInternalField(volVectorField& field) const
    {
        forAll(field, cellI)
        {
            const point& C = field.mesh().C()[cellI];
            field[cellI] = profile.velocityAt(C);
        }
    }

    template<template<class> class PatchField, class GeoMesh>
    void applyToBoundary
    (
        const word name, GeometricField<vector, PatchField, GeoMesh>& field
    ) const
    {
        const label boundaryI = field.mesh().boundaryMesh().findPatchID(name);
        forAll(field.boundaryField()[boundaryI], cellI)
        {
            const point& face = field.mesh().Cf().boundaryField()[boundaryI][cellI];
            field.boundaryField()[boundaryI][cellI] = profile.velocityAt(face);
        }
    }
};

#endif
