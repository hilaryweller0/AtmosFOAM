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

    const scalar deformationScale;
    
    vector streamfunctionAt
    (
        const label ip,
        const sphericalMeshData&,
        scalar time
    ) const;
};
}
#endif
