#ifndef DIVERGENT_VELOCITY_FIELD_H
#define DIVERGENT_VELOCITY_FIELD_H

#include "velocityField.H"

namespace Foam
{

class divergentVelocityField : public velocityField
{
    public:
    const bool applyProjection;
    void project(surfaceScalarField& phi) const;
    divergentVelocityField
    (
        const dictionary& dict,
        const bool applyProjection__ = false
    );
    virtual ~divergentVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi, scalar time) const;
    void applyToBoundary
    (
        surfaceScalarField& phi,
        const label patchI,
        scalar time
    ) const;
    virtual vector velocityAt(const point& p, scalar time) const = 0;
    
    protected:
    const dimensionedScalar endTime_;
};

}

#endif


