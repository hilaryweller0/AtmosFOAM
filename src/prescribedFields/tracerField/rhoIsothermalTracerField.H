#ifndef RHO_ISOTHERMAL_TRACER_FIELD_H
#define RHO_ISOTHERMAL_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class rhoIsothermalTracerField final : public tracerField
{
    public:
    TypeName("rhoIsothermal");
    rhoIsothermalTracerField
    (
        const dictionary& dict,
        const advectable& velocityField
    );

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar a_, rho0_, H_;
};

#endif


