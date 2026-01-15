#include "geodesicConeTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(geodesicConeTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicConeTracerField, dict);

geodesicConeTracerField::geodesicConeTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
R(readScalar(dict.lookup("radius"))*M_PI/180),
maxTracer(readScalar(dict.lookup("maxTracer"))),
lonCentre(readScalar(dict.lookup("lonCentre"))*M_PI/180),
latCentre(readScalar(dict.lookup("latCentre"))*M_PI/180)
{};

scalar geodesicConeTracerField::tracerAt(const point& p, const Time& t) const
{
    const polarPoint pp = convertToPolar(p);
    return maxTracer*(1-r(pp)/R);
}

scalar geodesicConeTracerField::r(const polarPoint& p) const
{
    return min(R, sqrt(sqr(p.lon()-lonCentre) + sqr(p.lat()-latCentre)));
}

}
