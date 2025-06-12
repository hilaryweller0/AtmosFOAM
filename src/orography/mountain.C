#include "mountain.H"

namespace Foam
{
defineTypeNameAndDebug(mountain, 0);
defineRunTimeSelectionTable(mountain, dict);

autoPtr<mountain> mountain::New(const dictionary& dict)
{
    const word mountainType(dict.lookup("type"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(mountainType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mountain::New(const dictionary&)"
        ) << "Unknown type "
          << mountainType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<mountain>
    (
       cstrIter()(dict)
    );
}

}
