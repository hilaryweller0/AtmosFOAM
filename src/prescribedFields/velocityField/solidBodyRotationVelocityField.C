#include "solidBodyRotationVelocityField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(solidBodyRotationVelocityField, 0);
addToRunTimeSelectionTable(velocityField, solidBodyRotationVelocityField, dict);

solidBodyRotationVelocityField::solidBodyRotationVelocityField(const dictionary& dict)
:
u0("speed", dimVelocity, dict.lookupOrDefault<scalar>("speed", scalar(10))),
z1("zeroVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("zeroVelocityHeight", scalar(4e3))),
z2("maxVelocityHeight", dimLength, dict.lookupOrDefault<scalar>("maxVelocityHeight", scalar(5e3)))
{};

vector solidBodyRotationVelocityField::streamfunctionAt
(
        const point& p,
        const Time& t
) const
{
    const vector unitNormal(0, -1, 0);
    const dimensionedScalar z("z", dimLength, p.z());
    const dimensionedScalar x("x", dimLength, p.x());
    const dimensionedScalar x1("x1", dimLength, scalar(-1000));
    const dimensionedScalar x2("x2", dimLength, scalar(1000));
    const dimensionedScalar w("angular_velocity", dimTime, scalar(100));

    dimensionedScalar psi("psi", cmptMultiply(dimVelocity, dimLength), scalar(0));
    //if (z.value() <= z1.value())
    //{
        // psi is zero
    //}
    //else if (z.value() <= z2.value())
    //{
        psi = 0.5*( z*z + x*x )/w;
    //}
    //else 
    //{
        // psi is zero
    //}

    return unitNormal * psi.value();
}

point solidBodyRotationVelocityField::initialPositionOf
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
        return point(p.x() - (u0*sqr(Foam::sin(0.5*M_PI*(z-z1)/(z2-z1)))*t).value(), p.y(), p.z());
    }
    else
    {
        return point(p.x() - u0.value()*t.value(), p.y(), p.z());
    }
}
