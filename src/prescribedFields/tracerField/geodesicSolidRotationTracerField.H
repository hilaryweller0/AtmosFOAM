#ifndef GEODESIC_SOLID_ROTATION_TRACER_FIELD_H
#define GEODESIC_SOLID_ROTATION_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"
#include "polarPoint.H"

namespace Foam
{
class geodesicSolidRotationTracerField final : public tracerField
{
    public:
    TypeName("geodesicSolidRotation");
    geodesicSolidRotationTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;

    // Private data
    const scalar magg;
    const scalar H0;
    const scalar Omega;
    const scalar u0;
};
}
#endif

