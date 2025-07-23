#include "deformationalNonDivergentVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(deformationalNonDivergentVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    deformationalNonDivergentVelocityField,
    dict
);

deformationalNonDivergentVelocityField::deformationalNonDivergentVelocityField
(
    const dictionary& dict
)
:
    nonDivergentVelocityField(dict),
    //deformationScale_(readScalar(dict.lookup("deformationScale"))),
    domainSize_(dict.lookup("domainSize")),
    backgroundFlow_(dict.lookupOrDefault<vector>("backgroundFlow", vector::zero))
{
    if (mag(domainSize_.z()) > SMALL)
    {
        FatalErrorIn("deformationalNonDivergentVelocityField")
            << " the size of the defomation domain non-zero in the x and y "
            << "direction (normal=k) but domainSize_ = " << domainSize_
            << exit(FatalError);
    }
};

vector deformationalNonDivergentVelocityField::streamfunctionAt
(
    const point& p,
    scalar time
) const
{
    const scalar T = endTime_.value();

    scalar xp = p.x() + domainSize_.x()/2 - backgroundFlow_[0]*time;
    scalar yp = p.y() + domainSize_.y()/2 - backgroundFlow_[1]*time;
    scalar C = backgroundFlow_[0]*domainSize_.x()/M_PI;
    
    return -vector(0,0,1)*
    (
        C*sqr
        (
            Foam::cos(M_PI*xp/domainSize_.x())
           *Foam::sin(M_PI*yp/domainSize_.y())
        )*Foam::cos(M_PI*time/T)
      + backgroundFlow_[0]*p.y()
      - backgroundFlow_[1]*p.x()
    );
}
}
