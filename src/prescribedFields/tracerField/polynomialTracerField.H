#ifndef POLYNOMIAL_TRACER_FIELD_H
#define POLYNOMIAL_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class polynomialTracerField final : public tracerField
{
    public:
    TypeName("polynomial");
    polynomialTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
};

#endif


