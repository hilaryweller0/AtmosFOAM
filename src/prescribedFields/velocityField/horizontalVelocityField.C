#include "horizontalVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(horizontalVelocityField, 0);
addToRunTimeSelectionTable(velocityField, horizontalVelocityField, dict);

horizontalVelocityField::horizontalVelocityField(const dictionary& dict)
:
    nonDivergentVelocityField(dict),
    u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
    z1
    (
        "zeroVelocityHeight",
        dimLength,
        dict.lookupOrDefault<scalar>("zeroVelocityHeight", scalar(4e3))
    ),
    z2
    (
        "maxVelocityHeight",
        dimLength,
        dict.lookupOrDefault<scalar>("maxVelocityHeight", scalar(5e3))
    )
{};

vector horizontalVelocityField::streamfunctionAt
(
        const point& p,
        scalar time
) const
{
    const vector unitNormal(0, -1, 0);
    const dimensionedScalar z("z", dimLength, p.z());

    dimensionedScalar psi("psi", cmptMultiply(dimVelocity, dimLength), scalar(0));
    if (z.value() <= z1.value())
    {
        // psi is zero
    }
    else if (z.value() <= z2.value())
    {
        psi = -0.5*u0*(z - z1 - (z2-z1)/M_PI*Foam::sin(M_PI*(z-z1)/(z2-z1)));
    }
    else 
    {
        psi = -0.5*u0*(2*z - z2 - z1);
    }

    return unitNormal * psi.value();
}

point horizontalVelocityField::initialPositionOf
(
    const point& p,
    scalar time
) const
{
    scalar z = p.z();

    if (z <= z1.value())
    {
        return p;
    }
    else if (z <= z2.value())
    {
        scalar ratio = (z - z1.value())/(z2.value() - z1.value());
        return point
        (
            p.x() - u0.value()*sqr(Foam::sin(0.5*M_PI*ratio))*time,
            p.y(),
            p.z()
        );
    }
    else
    {
        return point(p.x() - u0.value()*time, p.y(), p.z());
    }
}
}
