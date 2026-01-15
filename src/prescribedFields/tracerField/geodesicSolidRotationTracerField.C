#include "geodesicSolidRotationTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(geodesicSolidRotationTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicSolidRotationTracerField, dict);

geodesicSolidRotationTracerField::geodesicSolidRotationTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    magg(dimensionedScalar(dict.lookup("magg")).value()),
    H0(mag(dimensionedScalar(dict.lookup("H0")).value())),
    Omega(mag(dimensionedVector(dict.lookup("Omega")).value())),
    u0(mag(dimensionedScalar(dict.lookup("u0")).value()))
{};

scalar geodesicSolidRotationTracerField::tracerAt(const point& p, const Time& t) const
{
    const scalar a = mag(p);
    const scalar sinLat = p.z()/a;
    return H0 - (a*Omega*u0 + 0.5*sqr(u0))*sqr(sinLat)/magg;
}


}
