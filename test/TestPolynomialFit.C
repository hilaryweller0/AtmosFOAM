#include "TestPolynomialFit.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "TestableUpwindCorrFitData.H"

Test::PolynomialFit::PolynomialFit(
        const Foam::List<point>& stencilPoints,
        const Foam::label faceI)
{
    using namespace Foam;

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

    scalarList wts(stencilPoints.size(), scalar(1));
    const scalar unused_wLin = 0;

    extendedUpwindCellToFaceStencilNew* unusedStencil = NULL;
    bool linearCorrection = false;
    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;

    TestableUpwindCorrFitData
        fitData(mesh, *unusedStencil, linearCorrection, linearLimitFactor, centralWeight);
    fitData.calcFit(coefficients_, wts, stencilPoints, unused_wLin, faceI);
    coefficients_[0] += 1.0;
}

Test::PolynomialFit::PolynomialFit(const Foam::List<point>& stencilPoints)
{
    using namespace Foam;

    scalarList wts(stencilPoints.size(), scalar(1));
    const scalar unused_wLin = 0;
    const label unused_faceI = 0;

    fvMesh* unusedMesh = NULL;
    extendedUpwindCellToFaceStencilNew* unusedStencil = NULL;
    bool linearCorrection = false;
    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;
    const point p0(0, 0, 0);
    const bool pureUpwind = true;
    const vector idir(1, 0, 0);
    const vector jdir(0, 0, -1);
    const vector kdir(0, 1, 0);

/*    TestableUpwindCorrFitData
        fitData(*unusedMesh, *unusedStencil, linearCorrection, linearLimitFactor, centralWeight);
    fitData.calcFit(
            coefficients_,
            wts,
            stencilPoints,
            unused_wLin,
            p0,
            pureUpwind,
            idir,
            jdir,
            kdir,
            unused_faceI);*/
    coefficients_[0] += 1.0;
}

const Foam::scalarList Test::PolynomialFit::coefficients()
{
    return coefficients_;
}
