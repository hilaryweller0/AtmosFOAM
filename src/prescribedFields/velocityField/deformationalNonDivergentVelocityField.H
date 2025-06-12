#ifndef DEFORMATIONAL_NON_DIVERGENT_VELOCITY_FIELD_H
#define DEFORMATIONAL_NON_DIVERGENT_VELOCITY_FIELD_H

#include "nonDivergentVelocityField.H"

namespace Foam
{
class deformationalNonDivergentVelocityField final
     : public nonDivergentVelocityField
{
    public:
    TypeName("deformationalNonDivergent");
    deformationalNonDivergentVelocityField(const dictionary& dict);

    private:
    vector streamfunctionAt(const point& p, const Time& t) const;
    const scalar deformationScale_;
    const vector domainSize_;
};
}
#endif
