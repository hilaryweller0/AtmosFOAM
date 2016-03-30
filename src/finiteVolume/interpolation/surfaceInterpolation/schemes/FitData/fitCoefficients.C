#include "fitCoefficients.H"

Foam::fitCoefficients::fitCoefficients
(
    scalarList& coefficients,
    const List<point>& stencil,
    const bool linearCorrection,
    const scalar wLin
)
:
    coefficients(coefficients),
    linearCorrection(linearCorrection),
    wLin(wLin)
{
    coefficients.setSize(stencil.size());
}

scalar& Foam::fitCoefficients::operator[](int i)
{
    return coefficients[i];
}

scalar Foam::fitCoefficients::operator[](int i) const
{
    return coefficients[i];
}

label Foam::fitCoefficients::size()
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
