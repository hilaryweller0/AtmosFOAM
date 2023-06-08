#ifndef GEODESIC_NON_DIVERGENT_VELOCITY_FIELD_H
#define GEODESIC_NON_DIVERGENT_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"
#include "sphericalMeshData.H"

namespace Foam
{

class geodesicNonDivergentVelocityField : public velocityField
{
    public:
    geodesicNonDivergentVelocityField(const dictionary& dict);
    virtual ~geodesicNonDivergentVelocityField() {};
    scalar earthRadius() const { return earthRadius_; }

    private:
    scalar earthRadius_;

    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    scalar faceFlux(const face& f, const sphericalMeshData&, const Time& t) const;
    virtual vector streamfunctionAt
    (
        const label ip,
        const sphericalMeshData&,
        const Time& t
    ) const = 0;
    void project(surfaceScalarField& phi) const;
};

}

#endif


