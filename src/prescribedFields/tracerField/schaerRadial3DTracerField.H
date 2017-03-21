#ifndef SCHAER_RADIAL_3D_TRACER_FIELD_H
#define SCHAER_RADIAL_3D_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class schaerRadial3DTracerField final : public tracerField
{
    public:
    TypeName("schaerRadial3D");
    schaerRadial3DTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const point p0;
    const vector A;
};

#endif


