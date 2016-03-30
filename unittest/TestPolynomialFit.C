#include "TestPolynomialFit.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "TestableUpwindCorrFitData.H"
#include "AdaptivePolynomial.H"
#include "fitWeights.H"

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

    const scalar unused_wLin = 0;

    extendedUpwindCellToFaceStencilNew* unusedStencil = NULL;
    bool linearCorrection = false;
    scalar linearLimitFactor = 3.0;
    scalar centralWeight = 1e3;
    scalarList wts(stencilPoints.size(), scalar(1));

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

    scalar centralWeight = 1e3;
    const bool pureUpwind = false;

    fitWeights weights(stencilPoints.size());
    weights.setCentralWeight(centralWeight, pureUpwind);

    const scalar wLin = 0.6;

    scalar linearLimitFactor = 3.0;
    const point p0(0, 0, 0);
    const Basis basis(vector(1, 0, 0), vector(0, 1, 0), vector(0, 0, 1));
    const direction dimensions = 2;

    Foam::PolynomialFit<AdaptivePolynomial<cubicUpwindCPCFitPolynomial> > polynomialFit(
                linearCorrection,
                linearLimitFactor,
                centralWeight,
                dimensions);

    polynomialFit.fit(
            coefficients_,
            weights,
            stencilPoints,
            wLin,
            p0,
            pureUpwind,
            basis);

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
