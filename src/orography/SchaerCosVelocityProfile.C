#include "SchaerCosVelocityProfile.H"
#include "fvCFD.H"

SchaerCosVelocityProfile::SchaerCosVelocityProfile(const SchaerCosMountain& mountain, const IOdictionary& dict) :
    mountain(mountain),
    u0(readScalar(dict.lookup("maxVelocity"))),
    H(readScalar(dict.lookup("domainHeight")))
{};


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
    // TODO: worry about z over mountain
    scalar distanceToMountainStart = -mountain.a - p0.x();
    scalar timeToMountainStart = distanceToMountainStart / u0;
    scalar timeToCrossMountain = 4850.0; // TODO: calculate this from the big equation
    scalar timeAfterMountainEnd = t - timeToCrossMountain - timeToMountainStart;

    if (t <= timeToMountainStart) {
        // point not yet over mountain
        return point(p0.x() + u0 * t, 0, p0.z());
    } else if (t >= timeToMountainStart + timeToCrossMountain) {
        return point(p0.x() + distanceToMountainStart + 2.0*mountain.a + u0 * timeAfterMountainEnd, 0, p0.z());
    } else {
        // TODO: point is somewhere over the mountain, need to calculate this from the big equation
        return p0;
    }
}
