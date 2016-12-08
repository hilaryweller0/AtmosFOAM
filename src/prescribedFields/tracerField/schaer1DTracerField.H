#ifndef SCHAER_1D_TRACER_FIELD_H
#define SCHAER_1D_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class schaer1DTracerField final : public tracerField
{
    public:
    TypeName("schaer1D");
    schaer1DTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const point p0;
    const vector A;
};

#endif


