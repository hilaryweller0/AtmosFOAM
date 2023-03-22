#ifndef SPHERICAL_LAYER_TRACER_FIELD_H
#define SPHERICAL_LAYER_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class sphericalLayerTracerField final : public tracerField
{
    public:
    TypeName("sphericalLayer");
    sphericalLayerTracerField
    (
        const dictionary& dict,
        const advectable& velocityField
    );

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar a_, z1_, z2_, z0_, H_;
};

#endif


