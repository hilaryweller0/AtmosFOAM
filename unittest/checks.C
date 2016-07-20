#include "checks.H"

void check(Foam::fitCoefficients actual, Foam::fitCoefficients expected)
{
    forAll(actual, i)
    {
        CHECK(actual[i] == approx(expected[i]));
    }
}
