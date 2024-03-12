#ifndef SCHAER_RADIAL_TRACER_FIELD_H
#define SCHAER_RADIAL_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class schaerRadialTracerField final : public tracerField
{
    public:
    TypeName("schaerRadial");
    schaerRadialTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const point p0;
    const vector A;
    const label exponent;
};
}
#endif

