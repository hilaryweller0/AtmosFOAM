#ifndef BTF_VELOCITY_FIELD_H
#define BTF_VELOCITY_FIELD_H

#include "advectable.H"
#include "crossableMountain.H"
#include "nonDivergentVelocityField.H"
#include "btfTransform.H"

namespace Foam
{
enum class analyticSolution { HORIZONTAL_ONLY, VERTICAL_ONLY };

class btfVelocityField final : public nonDivergentVelocityField, public advectable
{
    public:
    TypeName("btf");
    btfVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    point initialPositionOf(const point& p, const Time& t) const;
    point initialHorizontalPositionOf(const point& p, const Time& t) const;

    const dimensionedScalar u0;
    const dimensionedScalar H;
    const autoPtr<mountain> m;
    const analyticSolution mode;
    const btfTransform transform; 
};
}
#endif

