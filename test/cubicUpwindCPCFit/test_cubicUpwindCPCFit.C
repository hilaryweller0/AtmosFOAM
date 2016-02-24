#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "TestableUpwindCorrFitData.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "cubicUpwindCPCFitPolynomial.H"

bool test()
{
    /*const fvMesh mesh
    (
        IOobject(
            "mesh",
            fileName("mesh"),
            Time("dummyTime", fileName("."), fileName("dummyCase")),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );*/

    const List<point> stencilPoints;
    scalarList coeffsi;
    scalarList wts(stencilPoints.size(), scalar(1));
    const scalar unused_wLin = 0;
    const label faceI = 0;

    const fvMesh* mesh = NULL;
    extendedUpwindCellToFaceStencilNew* unusedStencil = NULL;
    bool linearCorrection = false;
    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;

    TestableUpwindCorrFitData<cubicUpwindCPCFitPolynomial>
        fitData(*mesh, *unusedStencil, linearCorrection, linearLimitFactor, centralWeight);
    fitData.calcFit(coeffsi, wts, stencilPoints, unused_wLin, faceI);
    return false;
}

TEST_CASE("test")
{
    CHECK( test() );
}
