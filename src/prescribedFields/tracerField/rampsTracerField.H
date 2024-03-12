#ifndef RAMPS_TRACER_FIELD_H
#define RAMPS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class rampsTracerField final : public tracerField
{
    public:
    TypeName("ramps");
    rampsTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const label nRamps_;
    List<vector> centres_;
    List<vector> lengths_;
    scalarList maxTracers_;
};
}
#endif


