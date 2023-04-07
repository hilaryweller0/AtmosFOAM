#ifndef DIVERGENT_VELOCITY_FIELD_H
#define DIVERGENT_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"

namespace Foam
{

class divergentVelocityField : public velocityField
{
    public:
    const bool applyProjection;
    void project(surfaceScalarField& phi) const;
    divergentVelocityField(const bool applyProjection__ = false);
    virtual ~divergentVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    virtual vector velocityAt(const point& p, const Time& t) const = 0;
};

}

#endif


