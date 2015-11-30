#include "HorizontalVelocityProfile.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(HorizontalVelocityProfile, 0);
addToRunTimeSelectionTable(VelocityProfile, HorizontalVelocityProfile, dict);

HorizontalVelocityProfile::HorizontalVelocityProfile(const IOdictionary& dict) :
    u0(readScalar(dict.lookup("maxVelocity"))),
    z1(readScalar(dict.lookup("zeroVelocityHeight"))),
    z2(readScalar(dict.lookup("maxVelocityHeight")))
{};

HorizontalVelocityProfile::HorizontalVelocityProfile(
        const scalar u0,
        const scalar z1,
        const scalar z2) : 
    u0(u0),
    z1(z1),
    z2(z2)
{};

vector HorizontalVelocityProfile::velocityAt(const point& p) const
{
    const scalar z = p.z();

    if (z > z1 && z < z2)
    {
        return vector(u0*pow((Foam::sin(M_PI/2*(z-z1)/(z2-z1))),2), 0, 0);
    }
    else if (z >= z2)
    {
        return vector(u0, 0, 0);
    }
    else
    {
        return vector(0, 0, 0);
    }
}

scalar HorizontalVelocityProfile::streamFunctionAt(const point& p) const
{
    const scalar z = p.z();

    if (z <= z1) return 0;
    else if (z <= z2)
    {
        return -0.5*u0*(z - z1 - (z2-z1)/M_PI*Foam::sin(M_PI*(z-z1)/(z2-z1)));
    }
    else return -0.5*u0*(2*z - z2 - z1);
}

point HorizontalVelocityProfile::pointAtTime(const point& p0, const scalar t) const
{
    return point(p0.x() + u0 * t, 0, p0.z());
}
