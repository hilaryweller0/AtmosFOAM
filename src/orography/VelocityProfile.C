#include "VelocityProfile.H"

defineRunTimeSelectionTable(VelocityProfile, dict);

autoPtr<VelocityProfile> VelocityProfile::New(const IOdictionary& dict)
{
    const word velocityProfileType(dict.lookup("velocityProfileType"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(velocityProfileType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "VelocityProfile::New(const IOdictionary&)"
        ) << "Unknown velocityProfileType "
          << velocityProfileType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<VelocityProfile>
    (
       cstrIter()(dict)
    );
}
