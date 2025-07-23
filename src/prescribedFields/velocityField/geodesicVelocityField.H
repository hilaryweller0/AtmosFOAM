#ifndef GEODESIC_VELOCITY_FIELD_H
#define GEODESIC_VELOCITY_FIELD_H

#include "velocityField.H"

namespace Foam
{

class geodesicVelocityField : public velocityField
{
    protected:
    const dimensionedScalar earthRadius_;
    const dimensionedScalar endTime_;
        
    private:
    const bool applyProjection_;
    void applyToInternalField(surfaceScalarField& phi, scalar time) const;
    void applyToBoundary
    (
        surfaceScalarField& phi,
        const label patchI,
        scalar time
    ) const;
    virtual vector velocityAt(const point& p, scalar time) const = 0;

    public:
    void project(surfaceScalarField& phi) const;
    geodesicVelocityField(const dictionary& dict);
    virtual ~geodesicVelocityField() {};
    bool applyProjection() const { return applyProjection_; }
};

}

#endif


