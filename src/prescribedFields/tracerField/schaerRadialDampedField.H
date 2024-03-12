#ifndef SCHAER_RADIAL_DAMPED_FIELD_H
#define SCHAER_RADIAL_DAMPED_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class schaerRadialDampedField final : public tracerField
{
    public:
    TypeName("schaerDamped");
    schaerRadialDampedField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const string profileType;
    const scalar rho0;
    const scalar rhoAir;
    const point p0;
    const point p1;
    const vector A;
};
}
#endif

