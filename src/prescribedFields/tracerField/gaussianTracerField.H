#ifndef GAUSSIANS_TRACER_FIELD_H
#define GAUSSIANS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class gaussianTracerField final : public tracerField
{
    public:
    TypeName("gaussian");
    gaussianTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const dimensionedScalar dxMid;
    const point centre;
};

#endif


