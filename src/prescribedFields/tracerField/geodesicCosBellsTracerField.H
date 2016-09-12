#ifndef GEODESIC_COS_BELLS_TRACER_FIELD_H
#define GEODESIC_COS_BELLS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class geodesicCosBellsTracerField final : public tracerField
{
    public:
    TypeName("geodesicCosBells");
    geodesicCosBellsTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const dimensionedScalar R;
    const dimensionedScalar hmax;
    const dimensionedScalar b;
    const dimensionedScalar c;
};

#endif


