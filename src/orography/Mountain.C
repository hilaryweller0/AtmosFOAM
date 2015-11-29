#include "Mountain.H"
#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "addToRunTimeSelectionTable.H"

SchaerFineMountain::SchaerFineMountain(const scalar lambda) : lambda(lambda) {};

scalar SchaerFineMountain::heightAt(const point& p) const
{
        return sqr(Foam::cos(M_PI * p.x() / lambda));
}

defineRunTimeSelectionTable(Mountain, dict);

autoPtr<Mountain> Mountain::New(const IOdictionary& dict)
{
    const word mountainType(dict.lookup("mountainType"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(mountainType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Mountain::New(const IOdictionary&)"
        ) << "Unknown mountainType "
          << mountainType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<Mountain>
    (
       cstrIter()(dict)
    );
}
