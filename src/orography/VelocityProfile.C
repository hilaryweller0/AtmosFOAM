#include "VelocityProfile.H"
#include "HorizontalVelocityProfile.H"
#include "BtfVelocityProfile.H"
#include "Mountain.H"

VelocityProfile* VelocityProfile::lookup(IOdictionary dict)
{
    const word velocityProfileWord(dict.lookupOrDefault<word>("velocityProfile", "HORIZONTAL"));

    if (velocityProfileWord == "HORIZONTAL") {
        return new HorizontalVelocityProfile(dict);
    } else if (velocityProfileWord == "SCHAER_COS") {
        const SchaerCosMountain* mountain = new SchaerCosMountain(dict);
        return new BtfVelocityProfile(*mountain, dict);
    } else if (velocityProfileWord == "SCHAER_EXP") {
        const SchaerExpMountain* mountain = new SchaerExpMountain(dict);
        return new BtfVelocityProfile(*mountain, dict);
    } else {
        FatalErrorIn("VelocityProfile::lookup") << "Unknown velocityProfile '" << velocityProfileWord << "'" << exit(FatalError);
        return NULL;
    }
}
