#ifndef HORIZONTAL_VELOCITY_PROFILE_H
#define HORIZONTAL_VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "VelocityProfile.H"

class HorizontalVelocityProfile : public VelocityProfile
{
    public:
    TypeName("HorizontalVelocityProfile");
    HorizontalVelocityProfile(const IOdictionary& dict);
    HorizontalVelocityProfile(const scalar u0, const scalar z1, const scalar z2);
    vector velocityAt(const point& p) const;
    scalar streamFunctionAt(const point& p) const;
    point pointAtTime(const point& p0, const scalar t) const;
    private:
    const scalar u0;
    const scalar z1;
    const scalar z2;
};

#endif
