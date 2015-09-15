#include "BtfVelocityProfile.H"
#include "fvCFD.H"

BtfVelocityProfile::BtfVelocityProfile(const Mountain& mountain, const IOdictionary& dict) :
    mountain(mountain),
    u0(readScalar(dict.lookup("maxVelocity"))),
    H(readScalar(dict.lookup("domainHeight")))
{};

BtfVelocityProfile::BtfVelocityProfile(
        const Mountain& mountain,
        const scalar u0,
        const scalar H) : 
    mountain(mountain),
    u0(u0),
    H(H)
{};

vector BtfVelocityProfile::velocityAt(const point& p) const
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

scalar BtfVelocityProfile::streamFunctionAt(const point& p) const 
{
    scalar h = mountain.heightAt(p.x());
    return -u0 * H * (p.z() - h) / (H - h);
}

point BtfVelocityProfile::pointAtTime(const point& p0, const scalar t) const
{
    // TODO: worry about z over mountain
    scalar distanceToMountainStart = -mountain.halfWidth() - p0.x();
    scalar timeToMountainStart = distanceToMountainStart / u0;
    scalar timeAfterMountainEnd = t - mountain.timeToCross(u0, H) - timeToMountainStart;

    if (t <= timeToMountainStart) {
        // point not yet over mountain
        return point(p0.x() + u0 * t, 0, p0.z());
    } else if (t >= timeToMountainStart + mountain.timeToCross(u0, H)) {
        // point has passed over mountain
        return point(p0.x() + distanceToMountainStart + 2.0*mountain.halfWidth() + u0 * timeAfterMountainEnd, 0, p0.z());
    } else {
        // TODO: point is somewhere over the mountain, need to calculate this from the big equation
        return p0;
    }
}
