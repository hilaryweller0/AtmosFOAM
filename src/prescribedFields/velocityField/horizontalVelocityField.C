#include "horizontalVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(horizontalVelocityField, 0);
addToRunTimeSelectionTable(velocityField, horizontalVelocityField, dict);

horizontalVelocityField::horizontalVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
z1("zeroVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("zeroVelocityHeight", scalar(4e3))),
z2("maxVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("maxVelocityHeight", scalar(5e3)))
{};

vector horizontalVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
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
    const Time& t
) const
{
    const dimensionedScalar z("z", dimLength, p.z());

    if (z.value() <= z1.value())
    {
        return p;
    }
    else if (z.value() <= z2.value())
    {
        return point(p.x() - (u0*Foam::sin(0.5*M_PI*(z-z1)/(z2-z1))*t).value(), p.y(), p.z());
    }
    else
    {
        return point(p.x() - u0.value()*t.value(), p.y(), p.z());
    }
}
