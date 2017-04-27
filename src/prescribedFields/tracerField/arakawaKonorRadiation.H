#ifndef ARAKAWA_KONOR_RADIATION_TRACER_FIELD_H
#define ARAKAWA_KONOR_RADIATION_TRACER_FIELD_H

#include "tracerField.H"

class arakawaKonorRadiationTracerField final : public tracerField
{
    public:
    TypeName("arakawaKonorRadiation");
    arakawaKonorRadiationTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
};

#endif

