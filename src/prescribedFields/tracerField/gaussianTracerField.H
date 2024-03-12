#ifndef GAUSSIANS_TRACER_FIELD_H
#define GAUSSIANS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class gaussianTracerField final : public tracerField
{
    public:
    TypeName("gaussian");
    gaussianTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar width_;
    const point centre_;
    const scalar maxTracer_;
};
}
#endif


