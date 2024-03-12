#ifndef SCHAER_BLOCK_FIELD_H
#define SCHAER_BLOCK_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class schaerBlockField final : public tracerField
{
    public:
    TypeName("schaerBlock");
    schaerBlockField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar rho0;
    const scalar rhoAir;
    const point p0;
    const vector A;
};
}
#endif

