#include "polynomialTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(polynomialTracerField, 0);
addToRunTimeSelectionTable(tracerField, polynomialTracerField, dict);

polynomialTracerField::polynomialTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField)
{};

scalar polynomialTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    const scalar a_1 = 150;
    const scalar a_2 = 1e-4;
    const scalar a_3 = 1e-4;
    const scalar a_4 = 1e-4;
    const scalar a_5 = -1e-4;
    const scalar a_6 = 1e-4;
    const scalar a_7 = -1e-10;
    const scalar a_8 = 0;
    const scalar a_9 = 0;

    scalar x = p.x();
    scalar y = p.z();
    return a_1 + a_2 * x + a_3 * y + a_4 * sqr(x) + a_5 * x * y + a_6 * sqr(y) + a_7 * pow3(x) + a_8 * sqr(x) * y + a_9 * x * sqr(y);
}
