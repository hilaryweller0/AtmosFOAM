#ifndef THERMAL_PROFILE_H
#define THERMAL_PROFILE_H

#include "fvCFD.H"

class ThermalProfile
{
    public:
    ThermalProfile(
            const IOdictionary& environmentalProperties,
            const dimensionedVector g,
            const dimensionedScalar T0);
    scalar tracerAt(const point& p, const Time& t) const;
    vector gradAt(const point& p, const Time& t) const;

    private:
    const dimensionedVector g;
    const dimensionedScalar T0;
    const scalarList nLayers;
    const scalarList zLayers;
};

#endif
