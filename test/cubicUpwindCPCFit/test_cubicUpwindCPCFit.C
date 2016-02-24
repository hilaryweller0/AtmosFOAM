#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "catch.hpp"
#include "fvCFD.H"

#include "extendedUpwindCellToFaceStencilNew.H"
#include "TestableUpwindCorrFitData.H"
#include "PolynomialFit.H"

bool fitTo2DUniformMesh()
{
    const Time runTime(Time::controlDictName, fileName("."), fileName("."));
    const fvMesh mesh
    (
        IOobject(
            fvMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<point> stencilPoints(12, point(0, 0, 0));
    stencilPoints[0] = point(-0.5, 0, 0);
    stencilPoints[1] = point(0.5, 0, 0);
    stencilPoints[2] = point(-1.5, 0, 0);
    stencilPoints[3] = point(-2.5, 0, 0);
    stencilPoints[4] = point(-0.5, 0, 1);
    stencilPoints[5] = point(0.5, 0, 1);
    stencilPoints[6] = point(-1.5, 0, 1);
    stencilPoints[7] = point(-2.5, 0, 1);
    stencilPoints[8] = point(-0.5, 0, -1);
    stencilPoints[9] = point(0.5, 0, -1);
    stencilPoints[10] = point(-1.5, 0, -1);
    stencilPoints[11] = point(-2.5, 0, -1);
    scalarList coeffsi;
    scalarList wts(stencilPoints.size(), scalar(1));

    const scalar unused_wLin = 0;
    const label faceI = 11;

    extendedUpwindCellToFaceStencilNew* unusedStencil = NULL;
    bool linearCorrection = false;
    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;

    TestableUpwindCorrFitData
        fitData(mesh, *unusedStencil, linearCorrection, linearLimitFactor, centralWeight);
    fitData.calcFit(coeffsi, wts, stencilPoints, unused_wLin, faceI);
    coeffsi[0] += 1.0;
    //return coeffsi;
    return false;
}

TEST_CASE("fit full-size stencil to uniform 2D mesh")
{
    Test::PolynomialFit fit;
    CHECK(fit.coefficients()[0] == Approx(0.7));
}
