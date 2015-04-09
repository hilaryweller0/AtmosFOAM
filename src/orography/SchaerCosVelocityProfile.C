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
    scalar h = mountain.heightAt(p.x());
    return -u0 * H * (p.z() - h) / (H - h); // TODO: not sure why I need a negative in front here
}

point SchaerCosVelocityProfile::pointAtTime(const point& p0, const scalar t) const
{
    // TODO: worry about z over mountain
    scalar distanceToMountainStart = -mountain.a - p0.x();
    scalar timeToMountainStart = distanceToMountainStart / u0;
    scalar timeAfterMountainEnd = t - timeToCrossMountain() - timeToMountainStart;

    if (t <= timeToMountainStart) {
        // point not yet over mountain
        return point(p0.x() + u0 * t, 0, p0.z());
    } else if (t >= timeToMountainStart + timeToCrossMountain()) {
        // point has passed over mountain
        return point(p0.x() + distanceToMountainStart + 2.0*mountain.a + u0 * timeAfterMountainEnd, 0, p0.z());
    } else {
        // TODO: point is somewhere over the mountain, need to calculate this from the big equation
        return p0;
    }
}

scalar SchaerCosVelocityProfile::timeToCrossMountain() const
{
    const scalar x = 2.0 * mountain.a;
    const scalar alpha = M_PI/mountain.lambda;
    const scalar beta = M_PI/(2.0 * mountain.a);

    return x / u0 - mountain.h0/(4.0 * u0 * H) * 
        (x + 0.25*(Foam::sin(2.0 * (alpha + beta) * x) / (alpha + beta) +
         Foam::sin(2.0 * (alpha - beta) * x) / (alpha - beta)) +
         0.5 * (Foam::sin(2.0*alpha*x)/alpha + Foam::sin(2.0*beta*x)/beta));
}
