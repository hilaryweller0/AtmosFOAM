#ifndef VELOCITY_PROFILE_H
#define VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "Mountain.H"

class VelocityProfile
{
    public:
    virtual vector velocityAt(const point& p) const = 0;
    virtual ~VelocityProfile() {};
};

class SchaerCosVelocityProfile
{
    public:
    SchaerCosVelocityProfile(SchaerCosMountain& mountain, const scalar u0, const scalar H);
    vector velocityAt(const point& p) const;
    private:
    const SchaerCosMountain& mountain;
    const scalar u0;
    const scalar H;
};

#endif
