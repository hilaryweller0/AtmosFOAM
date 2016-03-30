#include "PolynomialFit.H"
#include "weightedMatrix.H"
#include "nullAdjuster.H"
#include "iterativeAdjuster.H"
#include "fitCoefficients.H"

template<class Polynomial>
Foam::PolynomialFit<Polynomial>::PolynomialFit
(
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight, 
    const direction dimensions
)
:
    linearCorrection_(linearCorrection),
    linearLimitFactor_(linearLimitFactor),
    centralWeight_(centralWeight),
    dim_(dimensions)
{}

template<class Polynomial>
autoPtr<Fit> Foam::PolynomialFit<Polynomial>::fit
(
    scalarList& coeffsi,
    fitWeights& weights,
    const List<point>& C,
    const scalar wLin,
    const point& origin,
    const bool pureUpwind,
    const Basis& basis
)
{
    fitCoefficients coefficients(coeffsi, C, linearCorrection_, wLin);

    const List<point> localStencil = toLocalCoordinates(C, origin, basis);
    Polynomial polynomial(localStencil, dim_);
    autoPtr<scalarRectangularMatrix> B = polynomial.matrix();
    weightedMatrix matrix(B());

    matrix.apply(weights);

    nullAdjuster adjuster
    (
        matrix, 
        coeffsi,
        weights
    );

/*    iterativeAdjuster adjuster
    (
        matrix, 
        coeffsi,
        weights,
        wLin,
        pureUpwind,
        linearCorrection_,
        linearLimitFactor_,
        centralWeight_
    );*/

    bool goodFit = adjuster.adjustWeights();
    coefficients.applyCorrection(goodFit);

    return autoPtr<Fit>(new Fit(
            C,
            coeffsi,
            goodFit,
            B->m()
    ));
}

template<class Polynomial>
List<point> Foam::PolynomialFit<Polynomial>::toLocalCoordinates
(
    const List<point>& stencilPoints,
    const point& origin,
    const Basis& basis)
{
    List<point> localPoints(stencilPoints.size(), point(0, 0, 0));

    scalar scale = scaleLocalCoordinates(origin, stencilPoints[0], basis);
    forAll(stencilPoints, i)
    {
        localPoints[i] = toLocalCoordinates(origin, stencilPoints[i], basis) / scale;
    }

    return localPoints;
}

template<class Polynomial>
point Foam::PolynomialFit<Polynomial>::toLocalCoordinates
(
    const point& origin,
    const point& p,
    const Basis& basis
)
{
    point d;

    d.x() = (p - origin)&basis.i;
    d.y() = (p - origin)&basis.j;
    #ifndef SPHERICAL_GEOMETRY
    d.z() = (p - origin)&basis.k;
    #else
    d.z() = mag(p) - mag(origin);
    #endif

    return d;
}

template<class Polynomial>
scalar Foam::PolynomialFit<Polynomial>::scaleLocalCoordinates
(
    const point& origin,
    const point& upwindPoint,
    const Basis& basis
)
{
    return cmptMax(cmptMag((toLocalCoordinates(origin, upwindPoint, basis))));
}
