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
    //
    //- The flux field from the streamfunction by application of Stokes' theorem
    //- for the internal flux field
    void applyToInternalField(surfaceScalarField& phi) const
    {
        const fvMesh& mesh = phi.mesh();
        forAll(phi, faceI)
        {
            // Circulate around the vertices of the face and sum to contribution
            // to the flux
            const face& f = mesh.faces()[faceI];
            point p0 = mesh.points()[f.last()];
            point p1 = mesh.points()[f.first()];
            point pmid = 0.5*(p0 + p1);
            phi[faceI] = profile.streamFunctionAt(pmid)*(p1.y() - p0.y());
            for(label ip = 1; ip < f.size(); ip++)
            {
                p0 = p1;
                p1 = mesh.points()[f[ip]];
                point pmid = 0.5*(p0 + p1);
                phi[faceI] += profile.streamFunctionAt(pmid)*(p1.y() - p0.y());
            }
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
    
    //- The flux field from the streamfunction by application of Stokes' theorem
    //- for the boundary
    void applyToBoundary(const word name, surfaceScalarField& phi) const
    {
        const fvMesh& mesh = phi.mesh();
        const label patchI = mesh.boundaryMesh().findPatchID(name);
        scalarField& bf = phi.boundaryField()[patchI];
        forAll(bf, faceI)
        {
            const face& f = mesh.boundaryMesh()[patchI][faceI];
            point p0 = mesh.points()[f.last()];
            point p1 = mesh.points()[f.first()];
            point pmid = 0.5*(p0 + p1);
            bf[faceI] = profile.streamFunctionAt(pmid)*(p1.y() - p0.y());
            for(label ip = 1; ip < f.size(); ip++)
            {
                p0 = p1;
                p1 = mesh.points()[f[ip]];
                point pmid = 0.5*(p0 + p1);
                bf[faceI] += profile.streamFunctionAt(pmid)*(p1.y() - p0.y());
            }
        }
    }
};

#endif
