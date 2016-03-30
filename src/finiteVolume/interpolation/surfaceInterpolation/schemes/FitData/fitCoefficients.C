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
