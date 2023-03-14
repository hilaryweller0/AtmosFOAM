#ifndef GEODESIC_SLOTTEDCYLINDER_TRACER_FIELD_H
#define GEODESIC_SLOTTEDCYLINDER_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class geodesicSlottedCylinderTracerField final : public tracerField
{
    public:
    TypeName("geodesicSlottedCylinder");
    geodesicSlottedCylinderTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar Rsphere;
    const scalar Rcylinder;
    const scalar hBackground;
    const scalar hmax;
};

#endif

