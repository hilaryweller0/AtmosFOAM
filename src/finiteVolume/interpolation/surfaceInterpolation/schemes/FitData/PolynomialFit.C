#include "PolynomialFit.H"
#include "weightedMatrix.H"
#include "nullAdjuster.H"
#include "iterativeAdjuster.H"

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
    const List<point> localStencil = toLocalCoordinates(C, origin, basis);
    Polynomial polynomial(localStencil, dim_);
    autoPtr<scalarRectangularMatrix> B = polynomial.matrix();
    weightedMatrix matrix(B());

    matrix.apply(weights);

    label stencilSize = C.size();
    coeffsi.setSize(stencilSize);

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
    applyCorrection(coeffsi, goodFit, wLin);

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

template<class Polynomial>
void Foam::PolynomialFit<Polynomial>::applyCorrection
(
    scalarList& coeffsi,
    const bool goodFit,
    const scalar wLin
)
{
    if (goodFit)
    {
        if (linearCorrection_)
        {
            // Remove the uncorrected linear coefficients
            coeffsi[0] -= wLin;
            coeffsi[1] -= 1 - wLin;
        }
        else
        {
            // Remove the uncorrected upwind coefficients
            coeffsi[0] -= 1.0;
        }
    }
    else
    {
        coeffsi = 0;
        
        if (linearCorrection_)
        {
            coeffsi[0] = 1-wLin;
            coeffsi[1] = -(1-wLin);
        }
    }
}
