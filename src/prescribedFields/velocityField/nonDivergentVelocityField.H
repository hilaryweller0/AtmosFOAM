#ifndef NON_DIVERGENT_VELOCITY_FIELD_H
#define NON_DIVERGENT_VELOCITY_FIELD_H

#include "velocityField.H"

namespace Foam
{

class nonDivergentVelocityField : public velocityField
{
    public:
    nonDivergentVelocityField(const dictionary& dict);
    virtual ~nonDivergentVelocityField() {};

    private:
    void applyToInternalField(surfaceScalarField& phi, scalar time) const;
    void applyToBoundary
    (
        surfaceScalarField& phi, const label patchI, scalar time
    ) const;
    scalar faceFlux(const face& f, const fvMesh& mesh, scalar time) const;
    virtual vector streamfunctionAt(const point& p, scalar time) const = 0;
    void project(surfaceScalarField& phi) const;
    
    protected:
    const dimensionedScalar endTime_;
};

}

#endif

