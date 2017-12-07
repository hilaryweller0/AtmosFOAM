#include "constructiveRotorsVelocityField.H"
#include "addToRunTimeSelectionTable.H"

#include "polarPoint.H"
#include "sphericalVector.H"

defineTypeNameAndDebug(constructiveRotorsVelocityField, 0);
addToRunTimeSelectionTable(velocityField, constructiveRotorsVelocityField, dict);

constructiveRotorsVelocityField::constructiveRotorsVelocityField(const dictionary& dict)
:
vmax(dict.lookupOrDefault<scalar>("maxVelocity", scalar(1))),
xmax(dict.lookupOrDefault<scalar>("xmax", scalar(1))),
ymax(dict.lookupOrDefault<scalar>("ymax", scalar(1))),
center(dict.lookupOrDefault<point>("centreOfSystem", point(0,0,0)))
{};

vector constructiveRotorsVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    scalar x = (p.x() - center.x())/xmax;
    scalar y = (p.y() - center.y())/ymax;
    scalar vx = vmax*Foam::sin(M_PI*x)*Foam::sin(M_PI*y);
    scalar vy = vmax*Foam::cos(1.5*M_PI*x)*Foam::cos(0.5*M_PI*y);
    return vector
    (
        vx,
        vy,
        0
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
