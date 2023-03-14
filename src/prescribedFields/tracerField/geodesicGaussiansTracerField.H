#ifndef GEODESIC_GAUSSIANS_TRACER_FIELD_H
#define GEODESIC_GAUSSIANS_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class geodesicGaussiansTracerField final : public tracerField
{
    public:
    TypeName("geodesicGaussians");
    geodesicGaussiansTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar R;
    const scalar hBackground;
    const scalar hmax;
    const scalar b;
};

#endif

