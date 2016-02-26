#ifndef Checks_H
#define Checks_H

#include "fvCFD.H"
#include "catch.hpp"

static Approx approx = Approx::custom().epsilon(0.001);

void check(Foam::scalarList actual, Foam::scalarList expected);

#endif

