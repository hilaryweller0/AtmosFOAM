#ifndef RADIAL_DISTANCE_H
#define RADIAL_DISTANCE_H

#include "advectable.H"
#include "tracerField.H"

class radialDistance final : public tracerField
{
    public:
    TypeName("radialDistance");
    radialDistance(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const point p0;
};

#endif


