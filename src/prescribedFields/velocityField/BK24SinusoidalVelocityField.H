#ifndef SINUSOIDAL_VELOCITY_FIELD_H
#define SINUSOIDAL_VELOCITY_FIELD_H

/*
A velocity field which is a cosine in each Cartesian direction. Returns:
    u = u0 + 0.5 u0 sin^2(pi xp/Lx) sin(2 pi yp/Ly) cos(pi t/T)
    v = v0 + 0.5 v0 sin^2(pi yp/Ly) sin(2 pi xp/Lx) cos(pi t/T)
    w = 0
    where xp = x + Lx/2 - u0 t and yp = y + Ly/2 - v0 t
*/

#include "divergentVelocityField.H"

namespace Foam
{
class BK24SinusoidalVelocityField final : public divergentVelocityField
{
    public:
    TypeName("BK24Sinusoidal");
    BK24SinusoidalVelocityField(const dictionary& dict);

    private:
    vector velocityAt(const point& p, scalar time) const;

    // The coeffiecients for different terms
    const scalar T_;
    const vector L_;
    const vector u0_;
};
}
#endif

