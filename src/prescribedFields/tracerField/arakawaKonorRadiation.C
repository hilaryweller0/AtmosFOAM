#include "arakawaKonorRadiation.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(arakawaKonorRadiationTracerField, 0);
addToRunTimeSelectionTable(tracerField, arakawaKonorRadiationTracerField, dict);

arakawaKonorRadiationTracerField::arakawaKonorRadiationTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField)
{};

scalar arakawaKonorRadiationTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    if (p.z() >= 3900 - SMALL && p.z() <= 4550 + SMALL)
    {
        return 10*Foam::sin(2*M_PI*p.x()/20e3);
    }
    else
    {
        return 0;
    }
}
