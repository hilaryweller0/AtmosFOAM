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
    if (mag(p.x()) <= 10e3)
    {
        if (p.z() >= 1300 - SMALL && p.z() <= 1950 + SMALL)
        {
            return -0.5*Foam::sin(2*M_PI*p.x()/20e3);
        }
        else if (p.z() >= 2600 - SMALL && p.z() <= 3250 + SMALL)
        {
            return 0.5*Foam::sin(2*M_PI*p.x()/20e3);
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
