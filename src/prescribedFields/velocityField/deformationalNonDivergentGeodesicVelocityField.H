#ifndef DEFORMATIONAL_NON_DIVERGENT_GEODESIC_VELOCITY_FIELD_H
#define DEFORMATIONAL_NON_DIVERGENT_GEODESIC_VELOCITY_FIELD_H

#include "geodesicNonDivergentVelocityField.H"

namespace Foam
{
class deformationalNonDivergentGeodesicVelocityField final : public geodesicNonDivergentVelocityField
{
    public:
    TypeName("deformationalNonDivergentGeodesic");
    deformationalNonDivergentGeodesicVelocityField(const dictionary& dict);

    private:

    const scalar radius;
    const scalar deformationScale;
    const scalar endTime;
    
    vector streamfunctionAt
    (
        const label ip,
        const sphericalMeshData&,
        const Time& t
    ) const;
};
}
#endif
