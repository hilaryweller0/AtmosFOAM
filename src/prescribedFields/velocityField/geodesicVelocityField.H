#ifndef GEODESIC_VELOCITY_FIELD_H
#define GEODESIC_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"

namespace Foam
{

class geodesicVelocityField : public velocityField
{
    // I need to create a mesh object with the spherical Cf and boundary Cf
    // But velocityField does not know about the mesh
    // Therefore the mesh object containing spherical data needs to be 
    // created by applyToInternalField and looked up when needed
    // Create the spherical data mesh object in another library (to replace 
    // spherical geometry)

    public:
    const bool applyProjection;
    void project(surfaceScalarField& phi) const;
    geodesicVelocityField(const bool applyProjection__ = false);
    virtual ~geodesicVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    virtual vector velocityAt(const point& p, const Time& t) const = 0;
};

}

#endif


