#ifndef DIVERGENT_VELOCITY_FIELD_H
#define DIVERGENT_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"

namespace Foam
{

class divergentVelocityField : public velocityField
{
    public:
    virtual ~divergentVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    virtual vector velocityAt(const point& p, const Time& t) const = 0;
};

}

#endif


