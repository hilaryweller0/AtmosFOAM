#include "constructiveRotorsVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"
#include "sphericalVector.H"

defineTypeNameAndDebug(constructiveRotorsVelocityField, 0);
addToRunTimeSelectionTable(velocityField, constructiveRotorsVelocityField, dict);

constructiveRotorsVelocityField::constructiveRotorsVelocityField(const dictionary& dict)
:
xmax(dict.lookupOrDefault<scalar>("xmax", scalar(1))),
ymax(dict.lookupOrDefault<scalar>("ymax", scalar(1))),
center(dict.lookupOrDefault<point>("centerOfSystem", point(0,0,0)))
{};

vector constructiveRotorsVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    scalar x = (p.x() - center.x())/xmax;
    scalar y = (p.z() - center.z())/ymax;
    return vector
    (
        Foam::sin(M_PI*x)*Foam::sin(M_PI*y),
        0,
        Foam::cos(1.5*M_PI*x)*Foam::cos(0.5*M_PI*y)
    );
}

point constructiveRotorsVelocityField::initialPositionOf
(
    const point& p,
    const Time& t
) const
{
    return p;
}
