#include "Checks.H"

void check(Foam::scalarList actual, Foam::scalarList expected)
{
    forAll(actual, i)
    {
        CHECK(actual[i] == approx(expected[i]));
    }
}
