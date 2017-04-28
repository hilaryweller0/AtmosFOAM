#ifndef ARAKAWA_KONOR_STRIPES_TRACER_FIELD_H
#define ARAKAWA_KONOR_STRIPES_TRACER_FIELD_H

#include "tracerField.H"

class arakawaKonorStripesTracerField final : public tracerField
{
    public:
    TypeName("arakawaKonorStripes");
    arakawaKonorStripesTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const scalar wavelength;
    const scalar z1Start;
    const scalar z1End;
    const scalar z2Start;
    const scalar z2End;
};

#endif

