#ifndef GEODESIC_VELOCITY_FIELD_H
#define GEODESIC_VELOCITY_FIELD_H

#include "velocityField.H"

namespace Foam
{

class geodesicVelocityField : public velocityField
{
    private:
    scalar earthRadius_;
    const bool applyProjection_;
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    virtual vector velocityAt(const point& p, const Time& t) const = 0;

    public:
    void project(surfaceScalarField& phi) const;
    geodesicVelocityField(const dictionary& dict);
    virtual ~geodesicVelocityField() {};
    scalar earthRadius() const { return earthRadius_; }
    bool applyProjection() const { return applyProjection_; }
};

}

#endif


