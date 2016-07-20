#ifndef checks_H
#define checks_H

#include "fvCFD.H"
#include "catch.hpp"
#include "fitCoefficients.H"

static Approx approx = Approx::custom().epsilon(0.001);

void check(Foam::fitCoefficients actual, Foam::fitCoefficients expected);

#endif
