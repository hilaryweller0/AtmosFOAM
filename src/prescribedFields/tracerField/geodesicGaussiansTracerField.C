#include "geodesicGaussiansTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(geodesicGaussiansTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicGaussiansTracerField, dict);

geodesicGaussiansTracerField::geodesicGaussiansTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    R(readScalar(dict.lookup("radius"))),
    hBackground(dict.lookupOrDefault<scalar>("hBackground", scalar(0))),
    hmax(readScalar(dict.lookup("hmax"))-hBackground),
    b(readScalar(dict.lookup("b")))
{};

scalar geodesicGaussiansTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar lon1 = 5.0*M_PI/6.0;
    scalar lat1 = 0;

    scalar lon2 = 7.0*M_PI/6.0;
	scalar lat2 = 0;

    const point centre1
    (
        R * Foam::cos(lat1) * Foam::cos(lon1),
        R * Foam::cos(lat1) * Foam::sin(lon1),
        R * Foam::sin(lat1)
    );

    const point centre2
    (
        R * Foam::cos(lat2) * Foam::cos(lon2),
        R * Foam::cos(lat2) * Foam::sin(lon2),
        R * Foam::sin(lat2)
    );

    return hBackground +
    (
        hmax * Foam::exp(-b*magSqr(p - centre1)/sqr(R))
      + hmax * Foam::exp(-b*magSqr(p - centre2)/sqr(R))
    );
}
