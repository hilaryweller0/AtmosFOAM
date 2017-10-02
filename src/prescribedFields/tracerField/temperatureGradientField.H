#ifndef TEMPERATURE_GRADIENT_FIELD_H
#define TEMPERATURE_GRADIENT_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class temperatureGradientField final : public tracerField
{
    public:
    TypeName("temperatureGradient");
    temperatureGradientField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar T0;
    const scalar g;
    const scalar cp;
    const scalar theta0;
    const point p0;
    const string gradientType;
};

#endif


