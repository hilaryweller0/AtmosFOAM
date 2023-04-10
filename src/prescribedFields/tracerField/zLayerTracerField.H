#ifndef Z_LAYER_TRACER_FIELD_H
#define Z_LAYER_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

class zLayerTracerField final : public tracerField
{
    public:
    TypeName("zLayer");
    zLayerTracerField
    (
        const dictionary& dict,
        const advectable& velocityField
    );

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar z1_, z2_, z0_, H_;
};

#endif


