/*
Creates a scalar field which is a cosine in each Cartesian direction. Intended
to force specific Fourier modes of a problem to become excited. Returns:

    f(x,y,z) = c + (a_x + b_x*cos(2pi*k_x*(x-phi_x)/L_x))
             * (a_y + b_y*cos(2pi*k_y*(y-phi_y)/L_y))
             * (a_z + b_z*cos(2pi*k_z*(z-phi_z)/L_z))
    Applied non zero in the box [xmin, xmax],[ymin, ymax],[zmin, zmax]
*/

#ifndef SINUSOIDAL_3D_TRACER_FIELD_H
#define SINUSOIDAL_3D_TRACER_FIELD_H

#include "advectable.H"
#include "tracerField.H"

namespace Foam
{
class sinusoidalTracerField final : public tracerField
{
    public:
    TypeName("sinusoidal");
    sinusoidalTracerField(const dictionary& dict, const advectable& velocityField);

    private:
    scalar tracerAt(const point& p, const Time& t) const;
    // The coeffiecients for different terms
    const scalar c, a_x, a_y, a_z,
                 b_x, b_y, b_z,
                 k_x, k_y, k_z,
                 L_x, L_y, L_z,
                 phi_x, phi_y, phi_z,
                 xmin, xmax, ymin, ymax, zmin, zmax;
    // The maximum and minimum possible values of the function
    const scalar funcMax, funcMin;
};
}
#endif

