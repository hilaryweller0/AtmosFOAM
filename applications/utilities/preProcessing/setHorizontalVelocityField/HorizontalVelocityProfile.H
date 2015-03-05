#ifndef HORIZONTAL_VELOCITY_PROFILE_H
#define HORIZONTAL_VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "VelocityProfile.H"

class HorizontalVelocityProfile : public VelocityProfile
{
    public:
    HorizontalVelocityProfile(const scalar u0, const scalar z1, const scalar z2);
    vector velocityAt(const point& p) const;
    scalar streamFunctionAt(const point& p) const;
    private:
    const scalar u0;
    const scalar z1;
    const scalar z2;
};

#endif
