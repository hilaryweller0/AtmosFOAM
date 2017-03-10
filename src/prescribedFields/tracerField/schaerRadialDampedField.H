#ifndef SCHAER_RADIAL_DAMPED_FIELD_H
#define SCHAER_RADIAL_DAMPED_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class schaerRadialDampedField final : public tracerField
{
    public:
    TypeName("schaerDamped");
    schaerRadialDampedField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const point p0;
    const vector A;
};

#endif


