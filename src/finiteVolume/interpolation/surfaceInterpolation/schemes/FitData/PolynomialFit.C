#include "PolynomialFit.H"
#include "FixedPolynomial.H"

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
    scalarList& wts,
    const List<point>& C,
    const scalar wLin,
    const point& origin,
    bool pureUpwind,
    const Basis& basis
)
{
    wts[0] = centralWeight_;
    if (!pureUpwind)
    {
        wts[1] = centralWeight_;
    }

    const List<point> localStencil = toLocalCoordinates(C, origin, basis);
    Polynomial polynomial(localStencil, dim_);
    WeightedMatrix matrix(polynomial.matrix());

    matrix.applyStencilPointWeights(wts);
    matrix.multiplyConstantAndLinearWeights(wts[0]);

    label stencilSize = C.size();
    coeffsi.setSize(stencilSize);

    bool goodFit = false;
    for (int iIt = 0; iIt < 8 && !goodFit; iIt++)
    {
        scalarRectangularMatrix Binv = matrix.pseudoInverse();

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[0]*wts[i]*Binv[0][i];
        }

        goodFit = isGoodFit(coeffsi, wLin);

        if (!goodFit)
        {
            increaseWeights(matrix, wts, pureUpwind, iIt == 0);
        }
    }

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
    return autoPtr<Fit>(new Fit(
            C,
            coeffsi,
            goodFit
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
    const Basis& basis)
{
    return cmptMax(cmptMag((toLocalCoordinates(origin, upwindPoint, basis))));
}

template<class Polynomial>
bool Foam::PolynomialFit<Polynomial>::isGoodFit
(
    const scalarList& coefficients, 
    const scalar wLin
)
{
    if (linearCorrection_)
    {
        return
            (mag(coefficients[0] - wLin) < linearLimitFactor_*wLin)
         && (mag(coefficients[1] - (1 - wLin)) < linearLimitFactor_*(1 - wLin))
         && eitherUpwindOrDownwindHasMaximumMagnitude(coefficients);
    }
    else
    {
        return (mag(coefficients[0] - 1.0) < linearLimitFactor_*1.0)
            && upwindCoefficientLargerThanSumOfOtherPositiveCoefficients(coefficients);
    }
}

template<class Polynomial>
bool Foam::PolynomialFit<Polynomial>::eitherUpwindOrDownwindHasMaximumMagnitude
(
    const scalarList& coefficients
)
{
    scalar maxCoeff = 0;
    label maxCoeffi = 0;
    for (label i=0; i<coefficients.size(); i++)
    {
        if (mag(coefficients[i]) > maxCoeff)
        {
            maxCoeff = mag(coefficients[i]);
            maxCoeffi = i;
        }
    }
    return maxCoeffi <= 1;
}

template<class Polynomial>
bool Foam::PolynomialFit<Polynomial>::upwindCoefficientLargerThanSumOfOtherPositiveCoefficients
(
    const scalarList& coefficients
)
{
    // go through all coeffsi except the first, add up all positive coeffs
    scalar positiveCoeffSum = 0;

    for(label i = 1; i < coefficients.size(); i++)
    {
        if (coefficients[i] > 0)
        {
            positiveCoeffSum += coefficients[i];
        }
    }

    return 3*coefficients[0] > positiveCoeffSum;
}

template<class Polynomial>
void Foam::PolynomialFit<Polynomial>::increaseWeights
(
    WeightedMatrix& matrix,
    scalarList& wts,
    bool pureUpwind,
    bool firstIteration
)
{
    wts[0] *= 10;
    if (linearCorrection_)
    {
        wts[1] *= 10;
    }
    else if (!pureUpwind && firstIteration)
    {
        wts[1] /= centralWeight_;
    }

    matrix.multiplyUpwindWeight(10);
    if (linearCorrection_)
    {
        matrix.multiplyDownwindWeight(10);
    }
    else if (!pureUpwind && firstIteration)
    {
        matrix.multiplyDownwindWeight(1.0/centralWeight_);
    }

    matrix.multiplyConstantAndLinearWeights(10);
}
