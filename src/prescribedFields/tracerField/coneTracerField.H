#ifndef CONES_TRACER_FIELD_H
#define CONES_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class coneTracerField final : public tracerField
{
    public:
    TypeName("cone");
    coneTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar width_;
    const point centre_;
    const scalar maxTracer_;
};

#endif


