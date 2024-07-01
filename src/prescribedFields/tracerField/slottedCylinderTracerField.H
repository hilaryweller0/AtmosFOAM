#ifndef SLOTTEDCYLINDER_TRACER_FIELD_H
#define SLOTTEDCYLINDER_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class slottedCylinderTracerField final : public tracerField
{
    public:
    TypeName("slottedCylinder");
    slottedCylinderTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar Rcylinder;
    const scalar Tbackground;
    const scalar Tmax;
    const point C1;
    const point C2;
    const scalar slotWidth;
};
}
#endif

