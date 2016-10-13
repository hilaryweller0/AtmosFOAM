#ifndef DEFORMATIONAL_NON_DIVERGENT_GEODESIC_VELOCITY_FIELD_H
#define DEFORMATIONAL_NON_DIVERGENT_GEODESIC_VELOCITY_FIELD_H

#include "nonDivergentVelocityField.H"

class deformationalNonDivergentGeodesicVelocityField final : public nonDivergentVelocityField
{
    public:
    TypeName("deformationalNonDivergentGeodesic");
    deformationalNonDivergentGeodesicVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    const dimensionedScalar radius;
};

#endif
