#ifndef GEODESIC_DEFORMATIONAL_DIVERGENT_VELOCITY_FIELD_H
#define GEODESIC_DEFORMATIONAL_DIVERGENT_VELOCITY_FIELD_H

#include "geodesicVelocityField.H"

class geodesicDeformationalDivergentVelocityField final : public geodesicVelocityField
{
    public:
    TypeName("geodesicDeformationalDivergent");
    geodesicDeformationalDivergentVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, const Time& t) const;
};

#endif


