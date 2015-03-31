#ifndef VELOCITY_PROFILE_H
#define VELOCITY_PROFILE_H

#include "fvCFD.H"

class VelocityProfile
{
    public:
    virtual vector velocityAt(const point& p) const = 0;
    virtual scalar streamFunctionAt(const point& p) const = 0;
    virtual point pointAtTime(const point& p0, const scalar t) const = 0;
    virtual ~VelocityProfile() {};
    static VelocityProfile* lookup(IOdictionary dict);
};

#endif
