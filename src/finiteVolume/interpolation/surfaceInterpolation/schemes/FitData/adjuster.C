#include "adjuster.H"

Foam::adjuster::adjuster
(
    polynomialMatrix& matrix, 
    scalarList& coefficients,
    scalarList& wts,
    const scalar wLin,
    const bool pureUpwind,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    matrix(matrix),
    coefficients(coefficients),
    wts(wts),
    wLin(wLin),
    pureUpwind(pureUpwind),
    linearCorrection(linearCorrection),
    linearLimitFactor(linearLimitFactor),
    centralWeight(centralWeight)
{}

bool Foam::adjuster::adjustWeights()
{
    bool goodFit = false;
    for (int iIt = 0; iIt < 8 && !goodFit; iIt++)
    {
        scalarRectangularMatrix Binv = matrix.pseudoInverse();

        for (label i=0; i<coefficients.size(); i++)
        {
            coefficients[i] = wts[0]*wts[i]*Binv[0][i];
        }

        goodFit = isGoodFit();

        if (!goodFit)
        {
            increaseWeights(iIt == 0);
        }
    }

    return goodFit;
}

bool Foam::adjuster::isGoodFit()
{
    if (linearCorrection)
    {
        return
            (mag(coefficients[0] - wLin) < linearLimitFactor*wLin)
         && (mag(coefficients[1] - (1 - wLin)) < linearLimitFactor*(1 - wLin))
         && eitherUpwindOrDownwindHasMaximumMagnitude();
    }
    else
    {
        return (mag(coefficients[0] - 1.0) < linearLimitFactor*1.0)
            && upwindCoefficientLargerThanSumOfOtherPositiveCoefficients();
    }
}

bool Foam::adjuster::eitherUpwindOrDownwindHasMaximumMagnitude()
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

bool Foam::adjuster::upwindCoefficientLargerThanSumOfOtherPositiveCoefficients()
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

void Foam::adjuster::increaseWeights(bool firstIteration)
{
    wts[0] *= 10;
    if (linearCorrection)
    {
        wts[1] *= 10;
    }
    else if (!pureUpwind && firstIteration)
    {
        wts[1] /= centralWeight;
    }

    matrix.multiplyUpwindWeight(10);
    if (linearCorrection)
    {
        matrix.multiplyDownwindWeight(10);
    }
    else if (!pureUpwind && firstIteration)
    {
        matrix.multiplyDownwindWeight(1.0/centralWeight);
    }

    matrix.multiplyConstantAndLinearWeights(10);
}
