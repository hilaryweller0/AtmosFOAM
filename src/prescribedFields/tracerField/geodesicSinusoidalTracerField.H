/*
Creates a scalar field which is a cosine in each spherical direction. Returns:
    f(x,y,z) = c + (a_x + b_x*Foam::cos(k_x*(x-phi_x)))
                 * (a_y + b_y*Foam::cos(k_y*(y-phi_y))
                 * (a_z + b_z*Foam::cos(2*M_PI*k_z*(z-phi_z)/L_z));
    Where x = longitude, y = latitude and z = altitude
*/

#ifndef GEODESIC_SINUSOIDAL_TRACER_FIELD_H
#define GEODESIC_SINUSOIDAL_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class geodesicSinusoidalTracerField final : public tracerField
{
    public:
    TypeName("geodesicSinusoidal");
    geodesicSinusoidalTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    const scalar Rsphere;
    const scalar c, a_x, a_y, a_z,
    b_x, b_y, b_z,
    k_x, k_y, k_z,
    L_z, phi_x, phi_y, phi_z;
};
}
#endif

