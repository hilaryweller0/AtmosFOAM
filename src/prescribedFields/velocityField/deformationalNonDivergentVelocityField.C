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
    deformationScale_(readScalar(dict.lookup("deformationScale"))),
    domainSize_(dict.lookup("domainSize"))
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
    const Time& time
) const
{
    const scalar T = time.endTime().value();
    const scalar t = time.value();

    return vector(0,0,1)*
    (
        deformationScale_*sqr(0.5*domainSize_.x()/M_PI)/T*sqr
        (
            Foam::sin(2*M_PI*(p.x()/domainSize_.x() - t/T))
           *Foam::cos(M_PI*p.y()/domainSize_.y())
        )
        *Foam::cos(M_PI*t/T)
      - domainSize_.x()*p.y()/T
    );
}
}
