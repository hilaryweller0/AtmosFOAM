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
hmax(dict.lookupOrDefault<scalar>("hmax", scalar(1))),
b(dict.lookupOrDefault<scalar>("b", scalar(0.1))),
c(dict.lookupOrDefault<scalar>("c", scalar(0.9)))
{};

scalar geodesicCosBellsTracerField::tracerAt(const point& p, const Time& t) const
{
    const polarPoint& polarp = convertToPolar(p);

    if (r1(polarp) < R/2)
    {
        return b + c*h1(polarp);
    }
    else if (r2(polarp) < R/2)
    {
        return b + c*h2(polarp);
    }
    else
    {
        return b;
    }
}

dimensionedScalar geodesicCosBellsTracerField::r1(const polarPoint& p) const
{
    return ri(p, 5*M_PI/6, 0);
}

dimensionedScalar geodesicCosBellsTracerField::r2(const polarPoint& p) const
{
    return ri(p, 7*M_PI/6, 0);
}

dimensionedScalar geodesicCosBellsTracerField::ri
(
        const polarPoint& p,
        const scalar centreLon,
        const scalar centreLat
) const
{
    return R * Foam::acos(Foam::sin(centreLat) * Foam::sin(p.lat()) + 
            Foam::cos(centreLat) * Foam::cos(p.lat()) * Foam::cos(p.lon() - centreLon));
}

scalar geodesicCosBellsTracerField::h1(const polarPoint& p) const
{
    return hi(p, r1(p));
}

scalar geodesicCosBellsTracerField::h2(const polarPoint& p) const
{
    return hi(p, r2(p));
}

scalar geodesicCosBellsTracerField::hi
(
        const polarPoint& p,
        const dimensionedScalar ri
) const
{
    return (hmax/2 * (1 + Foam::cos(M_PI * ri/(R/2)))).value();
}
