#include "slottedCylinderTracerField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(slottedCylinderTracerField, 0);
addToRunTimeSelectionTable(tracerField, slottedCylinderTracerField, dict);

slottedCylinderTracerField::slottedCylinderTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    Rcylinder(readScalar(dict.lookup("Rcylinder"))),
    Tbackground(dict.lookupOrDefault<scalar>("Tbackground", scalar(0))),
    Tmax(readScalar(dict.lookup("Tmax"))),
    C1(dict.lookup("C1")),
    C2(dict.lookup("C2")),
    slotWidth(readScalar(dict.lookup("slotWidth")))
{}

scalar slottedCylinderTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar r1 = mag(p - C1); 
    scalar r2 = mag(p - C2);
    scalar tracer = Tbackground;
    
    if
    (
        (r1 <= Rcylinder && (p.y() < C1.y() || mag(p.x() - C1.x()) > slotWidth))
     || (r2 <= Rcylinder && (p.y() < C2.y() || mag(p.x() - C2.x()) > slotWidth))
    )
    {
        tracer = Tmax;
    }
    
    return tracer;
}
}
