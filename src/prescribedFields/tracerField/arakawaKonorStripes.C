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
tracerField(velocityField)
{};

scalar arakawaKonorStripesTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    if (mag(p.x()) <= 100e3)
    {
        if (p.z() >= 900 - SMALL && p.z() <= 1900 + SMALL)
        {
            return -0.5*Foam::sin(2*M_PI*p.x()/200e3);
        }
        else if (p.z() >= 1900 - SMALL && p.z() <= 2900 + SMALL)
        {
            return 0.5*Foam::sin(2*M_PI*p.x()/200e3);
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
