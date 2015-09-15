#ifndef SCHAER_COS_VELOCITY_PROFILE_H
#define SCHAER_COS_VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "VelocityProfile.H"
#include "Mountain.H"

class BtfVelocityProfile : public VelocityProfile
{
    public:
    BtfVelocityProfile(const Mountain& mountain, const IOdictionary& dict);
    BtfVelocityProfile(const Mountain& mountain, const scalar u0, const scalar H);
    vector velocityAt(const point& p) const;
    scalar streamFunctionAt(const point& p) const;
    point pointAtTime(const point& p0, const scalar t) const;
    private:
    const Mountain& mountain;
    const scalar u0;
    const scalar H;
};

#endif
