#ifndef EXPONENTIAL_WITH_HEGIGHT_TRACER_FIELD_H
#define EXPONENTIAL_WITH_HEGIGHT_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class exponentialWithHeightTracerField final : public tracerField
{
    public:
    TypeName("exponentialWithHeight");
    exponentialWithHeightTracerField
    (
        const dictionary& dict,
        const advectable& velocityField
    );

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar z0_, T0_, H_;
};

#endif


