#ifndef GEODESIC_NON_DIVERGENT_VELOCITY_FIELD_H
#define GEODESIC_NON_DIVERGENT_VELOCITY_FIELD_H

#include "velocityField.H"
#include "sphericalMeshData.H"

namespace Foam
{

class geodesicNonDivergentVelocityField : public velocityField
{
    public:
    geodesicNonDivergentVelocityField(const dictionary& dict);
    virtual ~geodesicNonDivergentVelocityField() {};

    void applyToInternalField(surfaceScalarField& phi, scalar time) const;
    void applyToBoundary
    (
        surfaceScalarField& phi,
        const label patchI,
        scalar time
    ) const;
    scalar faceFlux(const face& f, const sphericalMeshData&, scalar time) const;
    virtual vector streamfunctionAt
    (
        const label ip,
        const sphericalMeshData&,
        scalar time
    ) const = 0;
    void project(surfaceScalarField& phi) const;
    
    protected:
    const dimensionedScalar earthRadius_;
    const dimensionedScalar endTime_;
};

}

#endif


