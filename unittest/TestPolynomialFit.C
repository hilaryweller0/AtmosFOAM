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

    TestableUpwindCorrFitData fitData(
            mesh,
            *unusedStencil,
            linearCorrection,
            linearLimitFactor,
            centralWeight);

    fitData.calcFit(coefficients_, wts, stencilPoints, unused_wLin, faceI);
    coefficients_[0] += 1.0;
}

Test::PolynomialFit::PolynomialFit(
        const Foam::List<point>& stencilPoints,
        bool linearCorrection)
{
    using namespace Foam;

    scalarList wts(stencilPoints.size(), scalar(1));
    const scalar wLin = 0.6;
    const label unused_faceI = 0;

    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;
    const point p0(0, 0, 0);
    const bool pureUpwind = false;
    const Basis basis(vector(1, 0, 0), vector(0, 0, -1), vector(0, 1, 0));
    const direction dimensions = 2;

    Foam::PolynomialFit<cubicUpwindCPCFitPolynomial> polynomialFit(
                linearCorrection,
                linearLimitFactor,
                centralWeight,
                dimensions);
    polynomialFit.fit(
            coefficients_,
            wts,
            stencilPoints,
            wLin,
            p0,
            pureUpwind,
            basis,
            unused_faceI);

    // TODO: shouldn't really have logic in the test code
    if (linearCorrection)
    {
        coefficients_[0] += wLin;
        coefficients_[1] += 1 - wLin;
    }
    else
    {
        coefficients_[0] += 1.0;
    }
}

const Foam::scalarList Test::PolynomialFit::coefficients()
{
    return coefficients_;
}
