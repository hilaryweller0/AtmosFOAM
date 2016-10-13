#ifndef NON_DIVERGENT_SPHERICAL_VELOCITY_FIELD_H
#define NON_DIVERGENT_SPHERICAL_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"

namespace Foam
{

class nonDivergentSphericalVelocityField : public velocityField
{
    public:
    virtual ~nonDivergentSphericalVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    scalar faceFlux(const face& f, const fvMesh& mesh, const Time& t) const;
    virtual vector streamfunctionAt(const point& p, const Time& t) const = 0;
};

}

#endif

