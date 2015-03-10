#ifndef SCHAER_COS_VELOCITY_PROFILE_H
#define SCHAER_COS_VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "VelocityProfile.H"
#include "Mountain.H"

class SchaerCosVelocityProfile : public VelocityProfile
{
    public:
    SchaerCosVelocityProfile(const SchaerCosMountain& mountain, const scalar u0, const scalar H);
    vector velocityAt(const point& p) const;
    scalar streamFunctionAt(const point& p) const;
    point pointAtTime(const point& p0, const scalar t) const;
    private:
    const SchaerCosMountain& mountain;
    const scalar u0;
    const scalar H;
};

#endif
