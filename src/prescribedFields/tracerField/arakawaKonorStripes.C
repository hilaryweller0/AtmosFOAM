#include "arakawaKonorStripes.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(arakawaKonorStripesTracerField, 0);
addToRunTimeSelectionTable(tracerField, arakawaKonorStripesTracerField, dict);

arakawaKonorStripesTracerField::arakawaKonorStripesTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(0.5))),
xOffset(dict.lookupOrDefault<scalar>("xOffset", scalar(0))),
wavelength(dict.lookupOrDefault<scalar>("wavelength", scalar(100e3))),
z1Start(dict.lookupOrDefault<scalar>("lowerWaveStart", scalar(1000))),
z1End(dict.lookupOrDefault<scalar>("lowerWaveEnd", scalar(2000))),
z2Start(dict.lookupOrDefault<scalar>("upperWaveStart", scalar(2000))),
z2End(dict.lookupOrDefault<scalar>("upperWaveEnd", scalar(3000)))
{};

scalar arakawaKonorStripesTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    scalar x = p.x() - xOffset;
    if (mag(x) <= wavelength / 2)
    {
        if (p.z() >= z1Start - 1 && p.z() <= z1End - 1)
        {
            return -rho0*Foam::sin(2*M_PI*x/wavelength);
        }
        else if (p.z() >= z2Start - 1 && p.z() <= z2End - 1)
        {
            return rho0*Foam::sin(2*M_PI*x/wavelength);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}
