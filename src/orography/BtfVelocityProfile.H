#ifndef BTF_VELOCITY_PROFILE_H
#define BTF_VELOCITY_PROFILE_H

#include "fvCFD.H"
#include "VelocityProfile.H"
#include "Mountain.H"

class BtfVelocityProfile : public VelocityProfile
{
    public:
    TypeName("BtfVelocityProfile");
    BtfVelocityProfile(const IOdictionary& dict);
    vector velocityAt(const point& p) const;
    scalar streamFunctionAt(const point& p) const;
    point pointAtTime(const point& p0, const scalar t) const;
    private:
    const autoPtr<Mountain> mountain;
    const scalar u0;
    const scalar H;
};

#endif
