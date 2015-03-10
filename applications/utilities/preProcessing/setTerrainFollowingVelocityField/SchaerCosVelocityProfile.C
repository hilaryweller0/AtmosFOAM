#include "SchaerCosVelocityProfile.H"
#include "fvCFD.H"

SchaerCosVelocityProfile::SchaerCosVelocityProfile(
        const SchaerCosMountain& mountain,
        const scalar u0,
        const scalar H) : 
    mountain(mountain),
    u0(u0),
    H(H)
{};

vector SchaerCosVelocityProfile::velocityAt(const point& p) const
{
    scalar h = mountain.heightAt(p.x());
    scalar dhdx = mountain.gradientAt(p.x());
    scalar u = H / (H - h);

    scalar w = H * dhdx * (H - p.z()) / pow(H - h, 2);

    if (p.z() >= h)
    {
        return vector(u0*u, 0, u0*w);
    }
    else
    {
        return vector(0,0,0);
    }
}

scalar SchaerCosVelocityProfile::streamFunctionAt(const point& p) const 
{
    return 0; // TODO
}

point SchaerCosVelocityProfile::pointAtTime(const point& p0, const scalar t) const
{
    return p0; // TODO
}
