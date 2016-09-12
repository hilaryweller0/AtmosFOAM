#include "geodesicCosBellsTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(geodesicCosBellsTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicCosBellsTracerField, dict);

geodesicCosBellsTracerField::geodesicCosBellsTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
R("radius", dimLength, dict.lookupOrDefault<scalar>("radius", scalar(6.3712e6))),
hmax("hmax", dimless, dict.lookupOrDefault<scalar>("hmax", scalar(1))),
b("b", dimless, dict.lookupOrDefault<scalar>("b", scalar(0.1))),
c("c", dimless, dict.lookupOrDefault<scalar>("b", scalar(0.9)))
{};

scalar geodesicCosBellsTracerField::tracerAt(const point& p, const Time& t) const
{
    return 0;
}
