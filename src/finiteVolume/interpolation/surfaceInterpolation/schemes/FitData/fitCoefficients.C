#include "fitCoefficients.H"

Foam::fitCoefficients::fitCoefficients
(
    const fitCoefficients& c
)
:
    coefficients(0),
    linearCorrection(c.linearCorrection),
    wLin(c.wLin)
{
    coefficients.append(c.coefficients);
}

Foam::fitCoefficients::fitCoefficients
(
    const label stencilSize,
    const bool linearCorrection,
    const scalar wLin
)
:
    coefficients(stencilSize, scalar(0)),
    linearCorrection(linearCorrection),
    wLin(wLin)
{}

void Foam::fitCoefficients::copyFrom(const fitCoefficients& source)
{
    forAll(source.coefficients, i)
    {
        coefficients[i] = source.coefficients[i];
    }
}

void Foam::fitCoefficients::copyInto(scalarList& target)
{
    target.append(coefficients);
}

scalar& Foam::fitCoefficients::operator[](int i)
{
    return coefficients[i];
}

scalar Foam::fitCoefficients::operator[](int i) const
{
    return coefficients[i];
}

label Foam::fitCoefficients::size() const
{
    return coefficients.size();
}

void Foam::fitCoefficients::applyCorrection(const bool goodFit)
{
    if (goodFit)
    {
        if (linearCorrection)
        {
            // Remove the uncorrected linear coefficients
            coefficients[0] -= wLin;
            coefficients[1] -= 1 - wLin;
        }
        else
        {
            // Remove the uncorrected upwind coefficients
            coefficients[0] -= 1.0;
        }
    }
    else
    {
        coefficients = 0;
        
        if (linearCorrection)
        {
            coefficients[0] = 1-wLin;
            coefficients[1] = -(1-wLin);
        }
    }
}

bool Foam::fitCoefficients::stable() const
{
    scalar upwind = coefficients[0];
    scalar downwind = coefficients[1];

    return mag(downwind) < upwind && upwind <= 1 + downwind && upwind >= 0.5 && central_are_largest();
}

bool Foam::fitCoefficients::central_are_largest() const
{
    scalar smallest_central_coefficient = min(coefficients[0], coefficients[1]);
    for (int i=2; i < coefficients.size(); i++)
    {
        if (mag(coefficients[i]) > smallest_central_coefficient) return false;
    }
    return true;
}

Foam::Ostream& Foam::operator<<
(
    Ostream& stream,
    const fitCoefficients& c
)
{
    stream << c.coefficients;
    return stream;
}
