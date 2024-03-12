#ifndef cosSqr_TRACER_FIELD_H
#define cosSqr_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class cosSqrTracerField final : public tracerField
{
    public:
    TypeName("cosSqr");
    cosSqrTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar halfWidth_;
    const point centre_;
    const scalar maxTracer_;
    const scalar backgroundTracer_;
};
}
#endif


