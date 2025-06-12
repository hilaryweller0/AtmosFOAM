#ifndef GEODESIC_COS_BELLS_TRACER_FIELD_H
#define GEODESIC_COS_BELLS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"
#include "polarPoint.H"

namespace Foam
{
class geodesicCosBellsTracerField final : public tracerField
{
    public:
    TypeName("geodesicCosBells");
    geodesicCosBellsTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    dimensionedScalar r1(const polarPoint& p) const;
    dimensionedScalar r2(const polarPoint& p) const;
    dimensionedScalar ri(const polarPoint& p, const scalar centreLon, const scalar centreLat) const;
    scalar h1(const polarPoint& p) const;
    scalar h2(const polarPoint& p) const;
    scalar hi(const polarPoint& p, const dimensionedScalar ri) const;

    const dimensionedScalar R;
    const scalar hmax;
    const scalar b;
    const scalar c;
};
}
#endif

