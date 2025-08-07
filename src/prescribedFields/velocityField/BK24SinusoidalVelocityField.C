#include "BK24SinusoidalVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(BK24SinusoidalVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    BK24SinusoidalVelocityField,
    dict
);

BK24SinusoidalVelocityField::BK24SinusoidalVelocityField
(
    const dictionary& dict
)
:
    divergentVelocityField(dict),
    T_(readScalar(dict.lookup("endTime"))),
    L_(dict.lookup("domainSize")),
    u0_(dict.lookup("backgroundFlow"))
{};

vector BK24SinusoidalVelocityField::velocityAt
(
    const point& p,
    scalar time
) const
{
    scalar xp = p[0] + 0.5*L_[0] - u0_[0]*time;
    scalar yp = p[1] + 0.5*L_[1] - u0_[1]*time;
    
    scalar u = u0_[0]
             + 0.5*u0_[0]*sqr(Foam::sin(M_PI*xp/L_[0]))
               *Foam::sin(2*M_PI*yp/L_[1])*Foam::cos(M_PI*time/T_);

    scalar v = u0_[1]
             + 0.5*u0_[1]*sqr(Foam::sin(M_PI*yp/L_[1]))
               *Foam::sin(2*M_PI*xp/L_[0])*Foam::cos(M_PI*time/T_);
               
    return vector(u,v,0);
}

}
