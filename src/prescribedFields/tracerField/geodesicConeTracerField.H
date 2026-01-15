#ifndef GEODESIC_CONE_TRACER_FIELD_H
#define GEODESIC_CONE_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"
#include "polarPoint.H"

namespace Foam
{
class geodesicConeTracerField final : public tracerField
{
    public:
    TypeName("geodesicCone");
    geodesicConeTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    scalar r(const polarPoint& p) const;

    const scalar R; // radius in radians, read in degrees
    const scalar maxTracer;
    const scalar lonCentre; // in radians, read in degrees
    const scalar latCentre; // in radians, read in degrees
};
}
#endif

