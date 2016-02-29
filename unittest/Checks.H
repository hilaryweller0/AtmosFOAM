#ifndef Checks_H
#define Checks_H

#include "fvCFD.H"
#include "catch.hpp"

static Approx approx = Approx::custom().epsilon(0.001);

void check(Foam::scalarList actual, Foam::scalarList expected);

template <size_t rows, size_t cols>
void check(
        const Foam::scalarRectangularMatrix& actual, 
        const scalar (&expected)[rows][cols]);

#ifdef NoRepository                                                             
#   include "Checks.C"
#endif

#endif

