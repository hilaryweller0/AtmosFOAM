#ifndef NON_DIVERGENT_VELOCITY_FIELD_H
#define NON_DIVERGENT_VELOCITY_FIELD_H

#include "fvCFD.H"

#include "velocityField.H"

namespace Foam
{

class nonDivergentVelocityField : public velocityField
{
    public:
    virtual ~nonDivergentVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    scalar faceFlux(const face& f, const fvMesh& mesh, const Time& t) const;
    virtual vector streamfunctionAt(const point& p, const Time& t) const = 0;
    void project(surfaceScalarField& phi) const;
};

}

#endif

