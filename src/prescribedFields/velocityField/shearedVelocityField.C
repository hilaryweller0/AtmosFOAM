#include "shearedVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(shearedVelocityField, 0);
addToRunTimeSelectionTable(velocityField, shearedVelocityField, dict);

shearedVelocityField::shearedVelocityField(const dictionary& dict)
:
    nonDivergentVelocityField(dict),
    gradU(dict.lookup("gradU")),
    u0(dict.lookup("referenceVelocity")),
    p0(dict.lookup("referencePosition")),
    acceleration(dict.lookupOrDefault<scalar>("acceleration", scalar(0))),
    gradUhat(gradU),
    normal(gradU ^ u0),
    magGradU(mag(gradU)),
    magU0(mag(u0)),
    u0hat(u0)
{
    scalar magNormal = mag(normal);
    if (magNormal < SMALL)
    {
        FatalErrorIn("shearedVelocityField::shearedVelocityField")
            << " the reference velocity and the velocity gradient should not be aligned but referenceVelocity = " << u0 << " and gradU = " << gradU << " and u0 ^ gradU = " << normal << exit(FatalError);
    }
    
    normal /= magNormal;
    gradUhat /= mag(gradU);
    u0hat /= mag(u0);
}

vector shearedVelocityField::streamfunctionAt
(
    const point& p,
    scalar time
) const
{
    scalar dist = (p - p0) & gradUhat;

    return normal*(1 + acceleration*time)*
    (
        0.5*sqr(dist)*magGradU
      + magU0*dist
    );
}

point shearedVelocityField::initialPositionOf
(
    const point& p,
    scalar time
) const
{
    vector U = u0 + (gradU & (p - p0))*u0hat;
    vector dist = U*time*(1 + 0.5*acceleration*time);
    return point(p - dist);
}
}
