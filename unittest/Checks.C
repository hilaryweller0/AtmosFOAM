#include "Checks.H"

void check(Foam::scalarList actual, Foam::scalarList expected)
{
    forAll(actual, i)
    {
        CHECK(actual[i] == approx(expected[i]));
    }
}

template <size_t rows, size_t cols>
void check
(
    const Foam::scalarRectangularMatrix& actual, 
    const scalar (&expected)[rows][cols]
)
{
    CHECK( actual.n() == rows );
    CHECK( actual.m() == cols );

    for (int i=0; i<actual.n(); i++)
    {
        for (int j=0; j<actual.m(); j++)
        {
            CHECK( actual[i][j] == approx(expected[i][j]) );
        }
    }
}

void checkStable(Foam::scalarList coefficients)
{
    scalar upwind = coefficients[0];
    scalar downwind = coefficients[1];

    CHECK( mag(downwind) < upwind );
    CHECK( upwind <= 1 + downwind );
}
