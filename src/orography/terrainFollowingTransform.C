#include "terrainFollowingTransform.H"

defineRunTimeSelectionTable(terrainFollowingTransform, dict);

autoPtr<terrainFollowingTransform> terrainFollowingTransform::New(const dictionary& dict)
{
    const word terrainFollowingTransformType(dict.lookup("type"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(terrainFollowingTransformType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "terrainFollowingTransform::New(const dictionary&)"
        ) << "Unknown type "
          << terrainFollowingTransformType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<terrainFollowingTransform>
    (
       cstrIter()(dict)
    );
}

